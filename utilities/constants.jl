# -----------------------------------------------------------------------------
# DEFINING CONSTANTS
# -----------------------------------------------------------------------------

# Gravitational Constants
  const μs = Dict(
  "Sun" => 1.32712440018e11,
  "Mercury" => 2.2032e4,
  "Venus" => 3.257e5,
  "Earth" => 3.986004415e5,
  "Moon" => 4.902799e3,
  "Mars" => 4.305e4,
  "Jupiter" => 1.266865361e8,
  "Saturn" => 3.794e7,
  "Uranus" => 5.794e6,
  "Neptune" => 6.809e6,
  "Pluto" => 9e2)

# Radii
  const rs = Dict(
  "Sun" => 696000.,
  "Mercury" => 2439.,
  "Venus" => 6052.,
  "Earth" => 6378.1363,
  "Moon" => 1738.,
  "Mars" => 3397.2,
  "Jupiter" => 71492.,
  "Saturn" => 60268.,
  "Uranus" => 25559.,
  "Neptune" => 24764.,
  "Pluto" => 1151.)

# Semi Major Axes
  const as = Dict(
  "Mercury" => 57909083.,
  "Venus" => 108208601.,
  "Earth" => 149598023.,
  "Moon" => 384400.,
  "Mars" => 227939186.,
  "Jupiter" => 778298361.,
  "Saturn" => 1429394133.,
  "Uranus" => 2875038615.,
  "Neptune" => 4504449769.,
  "Pluto" => 5915799000.)
