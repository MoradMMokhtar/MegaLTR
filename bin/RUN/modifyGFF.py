import re
import sys
import shutil

if len(sys.argv) != 4:
    print("Usage: python script.py input_file mapping_file order")
    sys.exit(1)

input_file = sys.argv[1]
mapping_file = sys.argv[2]
order = sys.argv[3]

mappings = {}

with open(mapping_file, 'r') as file:
    for line in file:
        word, replacement = line.strip().split('=')
        
        if order == '2':
            mappings[replacement.replace('>','')] = word.replace('>','')
        elif order == '1':
            mappings[word.replace('>','')] = replacement.replace('>','')
        else:
            print("order must be 1 or 2")
            sys.exit(1)
# print(mappings)       

temp_output_file = input_file + ".tmp"

with open(input_file, 'r') as infile, open(temp_output_file, 'w') as outfile:
    for line in infile:
        if not line.startswith('#'):
            for replacement, word in mappings.items():
                line = re.sub(re.escape(replacement) , word, line)
                # line = re.sub(re.escape(replacement) + r'\t', word+r'\t', line) ## use this if you have tab at the end of the ids

        outfile.write(line)


shutil.move(temp_output_file, input_file)

# print("Replacement complete!")
