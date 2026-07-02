use std::fmt;

#[derive(Debug)]
pub enum TrimError {
    SetpointsError(String),
    ParamsError(String),
}

impl fmt::Display for TrimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            TrimError::SetpointsError(msg) => write!(f, "invalid setpoints: {msg}"),
            TrimError::ParamsError(msg) => write!(f, "invalid params: {msg}"),
        }
    }
}

impl std::error::Error for TrimError {}
