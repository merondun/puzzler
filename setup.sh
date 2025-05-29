#!/bin/bash

echo "Adding Puzzler Pipeline to PATH..."

# Update the PATH variable with the correct directory
echo 'export PATH="$HOME/puzzler/bin:$PATH"' >> ~/.bashrc
chmod +x $HOME/puzzler/bin/*

# Apply changes immediately
source ~/.bashrc

echo "Installation complete! You can now run:"
echo "  puzzler_asm --sample --map_file"
echo "  puzzler_post --sample --map_file"
