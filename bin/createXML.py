import xml.etree.ElementTree as ET

def create_xml(input_file, template_file):
    # Parse the input and template XML files
    input_tree = ET.parse(input_file)
    template_tree = ET.parse(template_file)

    # Find the root elements in both trees
    input_root = input_tree.getroot()
    template_root = template_tree.getroot()

    # Remove existing MOVERS elements from the template
    for movers_element in template_root.findall('.//MOVERS'):
        template_root.remove(movers_element)

    # Iterate over the elements in the input root
    for element in input_root:
        # Append each element to the template root
        template_root.append(element)

    # Write the output to a file with indentation
    output_file = 'output.xml'
    ET.ElementTree(template_root).write(output_file, xml_declaration=True, encoding='utf-8', method="xml")

    print(f"XML file created successfully. Output file: {output_file}")

# Example usage
input_file = 'bin/input.xml'
template_file = 'bin/template.xml'
create_xml(input_file, template_file)
