use pyo3::{PyErr, exceptions::PyValueError};

/// error type for mavdac crate
#[derive(Debug, thiserror::Error)]
pub enum MavDACError {
    /// bad search pattern for images
    #[error("bad input pattern: {0}")]
    BadPattern(glob::PatternError),
    /// undreadable path (from glob)
    #[error("unreadable path: {0}")]
    UnreadablePath(glob::GlobError),
    /// io error wrapper
    #[error("{0}")]
    IOError(std::io::Error),
    /// fits image file is invalid
    #[error("Invalid FITs file: {0}")]
    InvalidFITS(String),
    /// invalid coordinate, e.g., out of bounds
    #[error("Coordinate error: {0}")]
    Coordinate(String),
    /// yaml file error
    #[error("YAML error: {0}")]
    YAMLError(serde_yaml::Error),
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
        PyValueError::new_err(format!("{value}"))
    }
}

impl From<serde_yaml::Error> for MavDACError {
    fn from(value: serde_yaml::Error) -> Self {
        MavDACError::YAMLError(value)
    }
}
