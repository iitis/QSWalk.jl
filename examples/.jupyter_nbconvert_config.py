## conversion command
# jupyter nbconvert --config .jupyter_nbconvert_config.py --to python *ipynb

## Extension of the file that should be written to disk
c.Exporter.file_extension = '.jl'

## This allows you to exclude input prompts from all templates if set to True.
c.TemplateExporter.exclude_input_prompt = True


