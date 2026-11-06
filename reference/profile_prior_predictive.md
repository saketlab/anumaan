# Prior Predictive Checks for a Joint Bernoulli Design

Draws `beta`/`eta` from the declared priors
(`design$stan_data$beta_sd`/`eta_sd` – i.e. whatever
[`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
was actually called with, not a separately-specified prior), computes
the implied joint-state distribution at each event's REAL covariates via
the SAME prediction kernel the Stan
model/[`profile_ppc`](https://saketlab.github.io/anumaan/reference/profile_ppc.md)/
[`profile_estimate`](https://saketlab.github.io/anumaan/reference/profile_estimate.md)
use, and summarizes it by group: per-class resistance margins, selected
pair probabilities, the two most extreme profiles (all-resistant,
all-susceptible), the resistance-count distribution, and a
concentration/entropy summary. This is a check on the PRIOR alone – it
never looks at `Y`, and a prior that looks implausible here is a prior
problem, not evidence the data disagrees with it.

## Usage

``` r
profile_prior_predictive(
  design,
  ids,
  n_draws = 1000L,
  ci_level = 0.95,
  seed = 123L
)
```

## Arguments

- design:

  List from
  [`profile_build_design`](https://saketlab.github.io/anumaan/reference/profile_build_design.md)
  (reads `$stan_data`: `X`, `B`, `V`, `beta_sd`, `eta_sd`; and
  `$design`: `class_names`).

- ids:

  Data frame, one row per event (matching `nrow(X)`), with a `group`
  column.

- n_draws:

  Integer. Prior draws. Default `1000L` (prior draws are cheap – no Stan
  sampling – so this can be large).

- ci_level:

  Numeric. Default `0.95`.

- seed:

  Integer. Default `123L`.

## Value

A list of data frames, each one row per group (x the family's own extra
key):

- `class_margins`:

  group x class: `mean`, `lower`, `upper` prior-predictive resistance
  probability.

- `pair_probabilities`:

  group x class-pair x cell (`RR`/`RS`/`SR`/`SS`): `mean`, `lower`,
  `upper`.

- `extreme_profiles`:

  group: `all_resistant_mean/ lower/upper`,
  `all_susceptible_mean/lower/upper`, and
  `independence_all_resistant`/`independence_all_susceptible` – what the
  all-R/all-S probability WOULD be under independence at this group's
  OWN mean marginal rates (the product of `class_margins`), so a
  prior-induced value far above the independence baseline signals strong
  prior-induced association, not just "a small number".

- `count_distribution`:

  group x `count` (`0..D` resistant classes): `mean`, `lower`, `upper`.

- `concentration`:

  group: `mean_entropy` (mean, across events and draws, of the per-event
  Shannon entropy of the `2^D`-state distribution, in nats),
  `max_entropy` (\\\log(2^D)\\, the uniform-distribution ceiling), and
  `entropy_ratio` (`mean_entropy / max_entropy` – near 0 means
  near-deterministic profiles; near 1 means close to the
  uniform/no-information prior).

## Details

Inspect whether the prior puts most of its mass on near-deterministic
profiles (entropy near 0) or implausibly strong cross-class associations
(e.g. all-resistant/all-susceptible probability far above what
independence at the marginal rates would give) across most hospitals
BEFORE fitting. When increasing `D` or the interaction order, re-run
this: many individually modest coefficients can combine into an extreme
combined score, so a prior that looked fine at a smaller `D` is not
guaranteed to still look fine at a larger one – rescaling by a fixed
formula without re-inspecting the induced distribution is exactly the
mistake this function exists to catch.
