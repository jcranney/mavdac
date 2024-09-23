use std::fmt;

use pyo3::{exceptions::PyValueError, PyErr};
pub type Result<T> = std::result::Result<T, MavDACError>;

#[derive(Debug)]
pub enum MavDACError {
    BadPattern(glob::PatternError),
    UnreadablePath(glob::GlobError),
    IOError(std::io::Error),
    InvalidFITS(String),
    LinalgError(String),
    Coordinate(String),
}

impl fmt::Display for MavDACError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MavDACError::BadPattern(..) => write!(f, "bad input pattern"),
            MavDACError::UnreadablePath(..) => write!(f, "unreadable path"),
            MavDACError::IOError(e) => write!(f, "{}", e.to_string()),
            MavDACError::InvalidFITS(s) => write!(f, "{}", s),
            MavDACError::LinalgError(e) => write!(f, "{}", e.to_string()),
            MavDACError::Coordinate(s) => write!(f, "{}", s),
        }
    }
}

impl std::error::Error for MavDACError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            MavDACError::BadPattern(err) => Some(err),
            MavDACError::UnreadablePath(err) => Some(err),
            MavDACError::IOError(err) => Some(err),
            MavDACError::InvalidFITS(_) => Some(self),
            MavDACError::LinalgError(_) => Some(self),
            MavDACError::Coordinate(_) => Some(self),
        }
    }
}

impl From<glob::PatternError> for MavDACError {
    fn from(value: glob::PatternError) -> Self {
        MavDACError::BadPattern(value)
    }
}

impl From<glob::GlobError> for MavDACError {
    fn from(value: glob::GlobError) -> Self {
        MavDACError::UnreadablePath(value)
    }
}

impl From<std::io::Error> for MavDACError {
    fn from(value: std::io::Error) -> Self {
        MavDACError::IOError(value)
    }
}

impl From<MavDACError> for PyErr {
    fn from(value: MavDACError) -> Self {
        match value {
            MavDACError::BadPattern(pattern_error) => PyValueError::new_err(pattern_error.to_string()),
            MavDACError::UnreadablePath(glob_error) => PyValueError::new_err(glob_error.to_string()),
            MavDACError::IOError(error) => PyValueError::new_err(error.to_string()),
            MavDACError::InvalidFITS(s) => PyValueError::new_err(s),
            MavDACError::LinalgError(s) => PyValueError::new_err(s),
            MavDACError::Coordinate(s) => PyValueError::new_err(s),
        }
    }
}