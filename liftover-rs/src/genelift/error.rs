use thiserror::Error;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Copy, Error)]
pub enum FeatureLiftError {
    #[error("error on sub-features")]
    SubFeatureError,
    #[error("no lift over candidates")]
    NoCandidates,
    #[error("too large lifted region: {0}")]
    TooLargeLiftedRegion(u64),
    #[error("too small lifted region: {0}")]
    TooSmallLiftedRegion(u64),
    #[error("multi-map")]
    MultiMap,
    #[error("multi-map sub feature")]
    MultiMapSubFeature,
    #[error("no common chromosome")]
    NoCommonChromosome,
    #[error("wrong exon order")]
    WrongExonOrder,
    #[error("wrong strand")]
    WrongStrand,
}

impl FeatureLiftError {
    pub fn error_message(&self) -> &'static str {
        match self {
            FeatureLiftError::SubFeatureError => "sub_feature_error",
            FeatureLiftError::NoCandidates => "no_lift_over",
            FeatureLiftError::TooLargeLiftedRegion(_) => "too_large_lifted_region",
            FeatureLiftError::TooSmallLiftedRegion(_) => "too_small_lifted_region",
            FeatureLiftError::MultiMap => "multi_map",
            FeatureLiftError::MultiMapSubFeature => "multi_map_sub_feature",
            FeatureLiftError::NoCommonChromosome => "no_common_chromosome",
            FeatureLiftError::WrongExonOrder => "wrong_exon_order",
            FeatureLiftError::WrongStrand => "wrong_strand",
        }
    }
}
