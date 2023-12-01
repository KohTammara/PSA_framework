import xml.etree.ElementTree as ET

def create_xml(input, template):
    # Parse the input and template XML files
    input_tree = ET.parse(input)
    template_tree = ET.parse(template)

    # Get the root element of the template tree
    template_root = template_tree.getroot()

    # Find the 'ROSETTASCRIPTS' element in the template tree
    rosetta_scripts = template_root.find('ROSETTASCRIPTS')

    # Add the input tree to the 'ROSETTASCRIPTS' element
    rosetta_scripts.append(input_tree.getroot())

    # Write the output to a file
    output_file = 'output.xml'
    template_tree.write(output_file)

    print(f"XML file created successfully. Output file: {output_file}")

# Example usage
input_file = 'input.xml'
template_file = 'template.xml'
create_xml(input_file, template_file)
