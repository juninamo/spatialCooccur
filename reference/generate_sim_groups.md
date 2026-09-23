# Generate multi-sample simulated data with disease-group structure

Wraps \[generate_sim()\] to produce N patients per group, each with a
(possibly noised) group-specific \`close_ratio\`, and optionally several
images (fields of view) per patient. Returns a single tidy data.frame
with \`x\`, \`y\`, \`cell_type\`, \`sample_id\`, \`group\`, \`patient\`
columns — directly consumable by \[nhood_enrichment_per_sample()\] and
the other \`\*\_per_sample\` helpers.

## Usage

``` r
generate_sim_groups(
  n_samples_per_group = 3,
  group_close_ratio = list(disease = 0.8, control = 0.2),
  n_types = 10,
  max_loc = 800,
  n_cells = 1500,
  test_type = "circle",
  distance_param = 50,
  between_sample_noise = 0.05,
  n_images_per_patient = 1,
  within_patient_noise = 0.05,
  seed = 1234
)
```

## Arguments

- n_samples_per_group:

  Integer, number of patients to generate per group.

- group_close_ratio:

  Named list of base \`close_ratio\` values, one entry per group, e.g.
  \`list(disease = 0.8, control = 0.2)\`.

- n_types, max_loc, test_type, distance_param:

  Passed through to \[generate_sim()\].

- n_cells:

  Number of cells per image; either a single number or a named list with
  one entry per group (e.g. to simulate groups whose images differ in
  size). When a group-specific value is given and \`max_loc\` is a
  single number, \`max_loc\` is scaled by \`sqrt(n_cells /
  n_cells_of_first_group)\` so that cell density is kept constant.

- between_sample_noise:

  SD of Gaussian noise added to the per-patient \`close_ratio\` around
  the group baseline (clipped to \[0,1\]).

- n_images_per_patient:

  Number of images per patient.

- within_patient_noise:

  SD of Gaussian noise added to the per-image \`close_ratio\` around the
  patient value (only used when \`n_images_per_patient \> 1\`).

- seed:

  Random seed (controls both the noise and the per-sample seeds passed
  to \`generate_sim\`).

## Value

A data.frame. \`sample_id\` identifies an image; \`patient\` identifies
the patient (equal to \`sample_id\` when \`n_images_per_patient = 1\`).

## Examples

``` r
df <- generate_sim_groups(n_samples_per_group = 2, n_images_per_patient = 2,
                          group_close_ratio = list(case = 0.8, control = 0.2),
                          n_types = 4, n_cells = 150, max_loc = 250,
                          test_type = "distribute", distance_param = 8)
unique(df[, c("sample_id", "patient", "group")])
#>                           sample_id   patient   group
#> case_1_img1_cell1       case_1_img1    case_1    case
#> case_1_img2_cell1       case_1_img2    case_1    case
#> case_2_img1_cell1       case_2_img1    case_2    case
#> case_2_img2_cell1       case_2_img2    case_2    case
#> control_1_img1_cell1 control_1_img1 control_1 control
#> control_1_img2_cell1 control_1_img2 control_1 control
#> control_2_img1_cell1 control_2_img1 control_2 control
#> control_2_img2_cell1 control_2_img2 control_2 control
```
