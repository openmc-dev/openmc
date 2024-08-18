# Argument list
$Env:args = ""

# Check for event-based
if ($Env:EVENT) {
  $Env:args = $Env:args + " --event "
}

# Run regression and unit tests
python --cov=openmc -v $Env:args tests