# bekkIRF 0.0.0.9003

- Add `higher_moment_mode` for KIRF and SIRF. The new default, `"hafner"`,
  averages raw higher moments before normalizing by averaged second moments;
  `"pathwise"` retains normalization within each simulated path. Carry the
  selected mode through bootstrap calculations and report it in IRF settings.
- Include the first worker error when all bootstrap IRF calculations fail.
- Correct scalar BEKK dynamics by converting the linear coefficients `a`, `b`,
  and `g` to square-root diagonal factors for the matrix recursion. Reject
  non-finite or negative scalar dynamics coefficients.
- Check bootstrap fit validity, finite parameters and log likelihood, and the
  iteration limit before accepting a draw. Retain fit diagnostics and report
  excluded draws instead of marking every returned fit as converged.
- Match the BEKKs joint-sign indicator at zero in both R and C++, preserving
  asymmetric responses to shocks with a zero component.
- Add regression coverage for these corrections. All five response types and
  the existing higher-moment normalization options remain available. KIRF
  methodological validation remains open.
