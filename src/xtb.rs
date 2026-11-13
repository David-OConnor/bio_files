//! Native GFN-FF molecular dynamics through xTB.

use std::{
    env,
    fs::{self, File},
    io::{self, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::{AtomicBool, Ordering},
    thread,
    time::{Duration, Instant},
};

use crate::{AtomGeneric, Xyz, md_out::create_run_dir, xyz::load_xyz_trajectory};

#[derive(Clone, Debug)]
pub struct XtbInput {
    pub atoms: Vec<AtomGeneric>,
    /// Picoseconds.
    pub dt: f64,
    pub num_steps: usize,
    pub temperature: f64,
    pub thermostat: bool,
    /// Steps between trajectory frames.
    pub output_interval: usize,
    pub relaxation_iters: Option<usize>,
    pub charge: f64,
    pub alpb_solvent: Option<String>,
    pub parameter_file: Option<PathBuf>,
    /// Otherwise XTB_EXECUTABLE, then xtb on PATH.
    pub executable: Option<PathBuf>,
    /// None permits an unbounded MD run.
    pub timeout: Option<Duration>,
}

#[derive(Clone, Debug)]
pub struct XtbFrame {
    /// Picoseconds, starting at zero.
    pub time: f64,
    pub geometry: Xyz,
    pub energy_hartree: Option<f64>,
}

#[derive(Debug)]
pub struct XtbOutput {
    pub frames: Vec<XtbFrame>,
    pub directory: PathBuf,
    pub log_text: String,
}

impl XtbInput {
    pub fn validate(&self) -> io::Result<()> {
        if self.atoms.is_empty()
            || self.atoms.iter().any(|a| {
                !(1..=86).contains(&a.element.atomic_number())
                    || a.posit.to_arr().iter().any(|v| !v.is_finite())
            })
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "GFN-FF requires finite coordinates and elements H–Rn",
            ));
        }

        if !self.dt.is_finite()
            || self.dt <= 0.
            || self.dt > 0.0005 * (1. + 1e-6)
            || self.num_steps == 0
            || self.num_steps > i32::MAX as usize
            || self.output_interval == 0
            || self.output_interval > i32::MAX as usize
            || !self.temperature.is_finite()
            || self.temperature <= 0.
            || self
                .relaxation_iters
                .is_some_and(|n| n == 0 || n > i32::MAX as usize)
            || self.timeout.is_some_and(|t| t.is_zero())
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Invalid xTB MD settings; use a positive time step of at most 0.0005 ps",
            ));
        }

        if !self.charge.is_finite() || (self.charge - self.charge.round()).abs() > 1e-6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "xTB requires a finite integer charge",
            ));
        }

        if self.alpb_solvent.as_deref().is_some_and(|s| {
            s.is_empty()
                || s.starts_with('-')
                || s.chars().any(|c| c.is_whitespace() || c.is_control())
        }) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "ALPB solvent must be a single xTB solvent name",
            ));
        }

        Ok(())
    }

    pub fn save(&self, directory: &Path) -> io::Result<()> {
        self.validate()?;
        fs::create_dir_all(directory)?;
        Xyz {
            atoms: self.atoms.clone(),
            comment: "GFN-FF MD".into(),
        }
        .save(&directory.join("geometry.xyz"))?;

        let mut input = File::create(directory.join("md.inp"))?;
        // xTB truncates time / step to an integer. Stay above the requested step count.
        let time = (self.num_steps as f64 + 0.001) * self.dt;
        let step = self.dt * 1000.;
        let dump = step * self.output_interval as f64 * (1. + 1e-12);
        writeln!(
            input,
            "$md\n time={time:.16}\n step={step:.16}\n dump={dump:.16}\n temp={}\n nvt={}\n hmass=0\n shake=0\n velo=false\n restart=false\n$end",
            self.temperature, self.thermostat
        )?;

        if let Some(iters) = self.relaxation_iters {
            writeln!(input, "$opt\n maxcycle={iters}\n$end")?;
        }

        if let Some(path) = &self.parameter_file {
            fs::copy(path, directory.join(".param_gfnff.xtb"))?;
        }

        fs::write(
            directory.join("command.txt"),
            format!("xtb {}\n", self.args().join(" ")),
        )?;

        Ok(())
    }

    fn args(&self) -> Vec<String> {
        let mut args: Vec<String> = ["geometry.xyz", "--gfnff", "--input", "md.inp", "--chrg"]
            .into_iter()
            .map(str::to_owned)
            .collect();
        args.push(self.charge.to_string());
        args.push(
            if self.relaxation_iters.is_some() {
                "--omd"
            } else {
                "--md"
            }
            .into(),
        );
        if let Some(solvent) = &self.alpb_solvent {
            args.extend(["--alpb".into(), solvent.clone()]);
        }

        args
    }

    /// Run in a new `md_out/xtb_{timestamp}/` directory; retain inputs, logs, topology, and
    /// restart data.
    pub fn run(&self) -> io::Result<XtbOutput> {
        self.run_with_cancel(&AtomicBool::new(false))
    }

    pub fn run_with_cancel(&self, cancel: &AtomicBool) -> io::Result<XtbOutput> {
        self.validate()?;
        let executable = self
            .executable
            .clone()
            .or_else(|| env::var_os("XTB_EXECUTABLE").map(PathBuf::from))
            .unwrap_or_else(|| "xtb".into());
        let executable = if executable.components().count() > 1 {
            fs::canonicalize(executable)?
        } else {
            executable
        };

        let directory = create_run_dir("xtb")?;
        self.save(&directory)?;
        let directory = fs::canonicalize(directory)?;
        let mut command = Command::new(&executable);
        command
            .current_dir(&directory)
            .args(self.args())
            .env("XTBPATH", &directory)
            .env("XTBHOME", &directory)
            .env("OMP_NUM_THREADS", "1")
            .env("MKL_NUM_THREADS", "1")
            .stdin(Stdio::null())
            .stdout(File::create(directory.join("stdout.log"))?)
            .stderr(File::create(directory.join("stderr.log"))?);

        #[cfg(windows)]
        {
            use std::os::windows::process::CommandExt;
            command.creation_flags(0x08000000);
        }

        let mut child = command.spawn().map_err(|e| {
            io::Error::new(
                e.kind(),
                format!("Unable to run xTB ({}): {e}", executable.display()),
            )
        })?;
        let start = Instant::now();
        let status = loop {
            if cancel.load(Ordering::Relaxed) || self.timeout.is_some_and(|t| start.elapsed() >= t)
            {
                let _ = child.kill();
                let _ = child.wait();
                let cancelled = cancel.load(Ordering::Relaxed);

                return Err(io::Error::new(
                    if cancelled {
                        io::ErrorKind::Interrupted
                    } else {
                        io::ErrorKind::TimedOut
                    },
                    format!(
                        "xTB {}; logs: {}",
                        if cancelled { "cancelled" } else { "timed out" },
                        directory.display()
                    ),
                ));
            }

            match child.try_wait() {
                Ok(Some(status)) => break status,
                Ok(None) => thread::sleep(Duration::from_millis(50)),
                Err(error) => {
                    let _ = child.kill();
                    let _ = child.wait();

                    return Err(error);
                }
            }
        };

        let log_text = format!(
            "{}\n{}",
            fs::read_to_string(directory.join("stdout.log"))?,
            fs::read_to_string(directory.join("stderr.log"))?
        );
        let lower = log_text.to_ascii_lowercase();
        if !status.success()
            || !directory.join("xtbmdok").is_file()
            || !lower.contains("normal termination of xtb")
            || lower.contains("emergency exit")
        {
            let tail: String = log_text
                .chars()
                .rev()
                .take(2000)
                .collect::<String>()
                .chars()
                .rev()
                .collect();

            return Err(io::Error::other(format!(
                "xTB MD failed ({status}); logs: {}\n{tail}",
                directory.display()
            )));
        }

        let frames = read_trajectory(
            &directory.join("xtb.trj"),
            self.dt * self.output_interval as f64,
        )?;
        let expected_frames = self.num_steps.div_ceil(self.output_interval);
        if frames.len() != expected_frames
            || frames.iter().any(|f| {
                f.geometry.atoms.len() != self.atoms.len()
                    || f.geometry
                        .atoms
                        .iter()
                        .zip(&self.atoms)
                        .any(|(a, b)| a.element != b.element)
            })
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Incomplete or mismatched xTB trajectory; logs: {}",
                    directory.display()
                ),
            ));
        }

        Ok(XtbOutput {
            frames,
            directory,
            log_text,
        })
    }
}

/// Read xTB's XYZ trajectory with `velo=false`. Frame spacing is in ps.
pub fn read_trajectory(path: &Path, frame_interval: f64) -> io::Result<Vec<XtbFrame>> {
    if !frame_interval.is_finite() || frame_interval <= 0. {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Invalid xTB frame interval",
        ));
    }

    load_xyz_trajectory(path)?
        .into_iter()
        .enumerate()
        .map(|(i, geometry)| {
            if geometry
                .atoms
                .iter()
                .any(|a| a.posit.to_arr().iter().any(|v| !v.is_finite()))
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Non-finite xTB trajectory coordinates",
                ));
            }

            let words: Vec<_> = geometry.comment.split_whitespace().collect();
            let energy_hartree = words
                .windows(2)
                .find(|w| w[0].eq_ignore_ascii_case("energy:"))
                .and_then(|w| w[1].replace(['D', 'd'], "E").parse::<f64>().ok())
                .or_else(|| geometry.comment.trim().parse::<f64>().ok());
            if energy_hartree.is_some_and(|e| !e.is_finite()) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Non-finite xTB trajectory energy",
                ));
            }

            Ok(XtbFrame {
                time: i as f64 * frame_interval,
                geometry,
                energy_hartree,
            })
        })
        .collect()
}
