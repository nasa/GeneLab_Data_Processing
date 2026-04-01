input=$1

# Remove leading and trailing brackets
input=$(echo "$input" | sed 's/^\[//; s/\]$//')

# Convert the input to KEY="VALUE" format with sanitized keys
interim_output=$(echo "$input" | sed -E 's/:/=/g' | tr ',' '\n')

while IFS='=' read -r key value; do
    # Remove leading and ending quotes from the value, if present
    value=$(echo "$value" | sed 's/^"\(.*\)"$/\1/')

    # Wrap the value in double quotes
    value="'$value'"

    # Remove leading spaces
    key=$(echo "$key" | sed 's/^ *//g')

    # Sanitize the key to follow valid Bash variable name format
    # Convert upper case to lower case, replace spaces with underscores, and remove square brackets
    sanitized_key=$(echo "$key" | tr '[:upper:] ' '[:lower:]_' | sed 's/[][]//g')


    # Append the key-value pair to the formatted output
    formatted_output+="$sanitized_key=$value\n"
done <<< "$interim_output"

echo -e "$formatted_output"
