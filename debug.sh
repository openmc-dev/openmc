./build.sh
echo "OpenMC is now done building and updating."
echo "Running ATR test..."
cd examples/atr/
python3 simulation.py