import re
import random
import sys
import shutil

def generate_random_word(length):
    return ''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789') for _ in range(length))

if len(sys.argv) != 4:
    print("Usage: python script.py input_file regex_pattern mapping_file")
    sys.exit(1)

input_file = sys.argv[1]
regex_to_replace = sys.argv[2]
mapping_file = sys.argv[3]

temp_output_file = input_file + ".tmp"

with open(input_file, 'r') as infile:
    with open(mapping_file, 'w') as mapfile, open(temp_output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                matched_words = set(re.findall(regex_to_replace, line))

                mappings = {}

                for word in matched_words:
                    replacement = '>' + generate_random_word(10)
                    mappings[word] = replacement

                    line = re.sub(re.escape(word) + r'\b', replacement, line)

                for word, replacement in mappings.items():
                    mapfile.write(f"{word}={replacement}\n")

                outfile.write(line)
            else:

                outfile.write(line)

shutil.move(temp_output_file, input_file)

# print("Replacement complete!")
