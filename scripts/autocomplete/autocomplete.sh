#!/usr/bin/env bash

CLI_COMMAND="edms"
RC_FILE=""
SHELL_TYPE="$(basename "$SHELL")"

# Determine shell and config file
if [[ "$SHELL_TYPE" == "bash" ]]; then
    RC_FILE="$HOME/.bashrc"
    AUTOCOMP_LINE="eval \"\$(register-python-argcomplete $CLI_COMMAND)\""
elif [[ "$SHELL_TYPE" == "zsh" ]]; then
    RC_FILE="$HOME/.zshrc"
    AUTOCOMP_LINE="eval \"\$(register-python-argcomplete $CLI_COMMAND)\""
    COMPINIT_LINE="autoload -U compinit && compinit"
else
    echo "Unsupported shell: $SHELL_TYPE"
    exit 1
fi

# Add compinit if needed (for Zsh)
if [[ "$SHELL_TYPE" == "zsh" ]]; then
    if ! grep -q "compinit" "$RC_FILE"; then
        echo "Adding compinit to $RC_FILE..."
        echo -e "\n$COMPINIT_LINE" >> "$RC_FILE"
    fi
fi

# Add autocomplete line if not already present
if ! grep -Fq "$AUTOCOMP_LINE" "$RC_FILE"; then
    echo "Adding autocomplete line to $RC_FILE..."
    echo -e "\n$AUTOCOMP_LINE" >> "$RC_FILE"
else
    echo "Autocomplete already enabled in $RC_FILE."
fi

echo "Done. Please run: source $RC_FILE"