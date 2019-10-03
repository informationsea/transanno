use failure::Fail;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Copy, Fail)]
pub enum FeatureLiftError {
    #[fail(display = "error on sub-features")]
    SubFeatureError,
    #[fail(display = "no lift over candidates")]
    NoCandidates,
    #[fail(display = "too large lifted region: {}", _0)]
    TooLargeLiftedRegion(u64),
    #[fail(display = "too small lifted region: {}", _0)]
    TooSmallLiftedRegion(u64),
    #[fail(display = "multi-map")]
    MultiMap,
    #[fail(display = "multi-map sub feature")]
    MultiMapSubFeature,
    #[fail(display = "no common chromosome")]
    NoCommonChromosome,
    #[fail(display = "wrong exon order")]
    WrongExonOrder,
    #[fail(display = "wrong strand")]
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
