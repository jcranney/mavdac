use std::fmt;

use pyo3::{exceptions::PyValueError, PyErr};
pub type Result<T> = std::result::Result<T, MavDACError>;

/// error type for mavdac crate
#[derive(Debug)]
pub enum MavDACError {
    /// bad search pattern for images
    BadPattern(glob::PatternError),
    /// undreadable path (from glob)
    UnreadablePath(glob::GlobError),
    /// io error wrapper
    IOError(std::io::Error),
    /// fits image file is invalid
    InvalidFITS(String),
    /// invalid coordinate, e.g., out of bounds
    Coordinate(String),
    /// yaml file error
    YAMLError(serde_yaml::Error),
}

impl fmt::Display for MavDACError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MavDACError::BadPattern(..) => write!(f, "bad input pattern"),
            MavDACError::UnreadablePath(..) => write!(f, "unreadable path"),
            MavDACError::IOError(e) => write!(f, "{}", e),
            MavDACError::InvalidFITS(s) => write!(f, "{}", s),
            MavDACError::Coordinate(s) => write!(f, "{}", s),
            MavDACError::YAMLError(e) => write!(f, "{}", e),
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
            MavDACError::Coordinate(_) => Some(self),
            MavDACError::YAMLError(err) => Some(err),
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
            MavDACError::Coordinate(s) => PyValueError::new_err(s),
            MavDACError::YAMLError(error) => PyValueError::new_err(error.to_string()),
        }
    }
}

impl From<serde_yaml::Error> for MavDACError {
    fn from(value: serde_yaml::Error) -> Self {
        MavDACError::YAMLError(value)
    }
}