import re
import sys
import shutil

if len(sys.argv) != 3:
    print("Usage: python script.py input_file mapping_file")
    sys.exit(1)

input_file = sys.argv[1]
mapping_file = sys.argv[2]

mappings = {}

with open(mapping_file, 'r') as file:
    for line in file:
        word, replacement = line.strip().split('=')
        mappings[replacement] = word

temp_output_file = input_file + ".tmp"

with open(input_file, 'r') as infile, open(temp_output_file, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):
            for replacement, word in mappings.items():
                line = re.sub(re.escape(replacement) + r'\b', word, line)

        outfile.write(line)


shutil.move(temp_output_file, input_file)

print("Replacement complete!")
