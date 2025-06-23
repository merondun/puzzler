#!/bin/bash

echo "Adding Puzzler Pipeline to PATH..."

# Update the PATH variable with the correct directory
INSTALL_PATH=$(dirname "$(realpath bin/puzzler)")
grep -qxF "export PATH=\$PATH:$INSTALL_PATH" ~/.bashrc || echo "export PATH=\$PATH:$INSTALL_PATH" >> ~/.bashrc
export PATH="$PATH:$INSTALL_PATH"
chmod +x bin/puzzler

# Apply changes immediately
source ~/.bashrc

echo "Installation complete! You can now run:"
echo "  puzzler --sample $SAMPLE --map_file examples/samples.tsv"
