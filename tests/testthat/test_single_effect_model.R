# =============================================================================
# Single Effect Model tests
#
# The original tests compared mvsusieR's R6 SingleEffectModel against susieR's
# single_effect_regression at the SER level. Since susieR's data representation
# changed (individual/ss classes vs DenseData/SSData R6), SER-level comparison
# requires incompatible data object setups.
#
# These tests are superseded by:
#   - test_r1_susieR_identity.R: compares full IBSS algorithm for R=1
#   - test_reference_identity.R: compares S3 vs R6 for R>1
# =============================================================================
