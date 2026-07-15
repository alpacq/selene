//! Error types for the selene program.

use std::fmt;

#[derive(Debug)]
pub enum TrimError {
    SetpointsError(String),
    ParamsError(String),
    ConvergenceError(String),
}

impl fmt::Display for TrimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TrimError::SetpointsError(msg) => write!(f, "invalid setpoints: {msg}"),
            TrimError::ParamsError(msg) => write!(f, "invalid params: {msg}"),
            TrimError::ConvergenceError(msg) => write!(f, "convergence error: {msg}"),
        }
    }
}

impl std::error::Error for TrimError {}

#[derive(Debug)]
pub enum LinearizationError {
    ConvergenceError(String),
}

impl fmt::Display for LinearizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LinearizationError::ConvergenceError(msg) => write!(f, "convergence error: {msg}"),
        }
    }
}

impl std::error::Error for LinearizationError {}
