from jinja2 import Template
from sys import argv

types = [
    {'name': 'rpe', 'code': 'type(rpe_var)'},
    {'name': 'single', 'code': 'real(sp)'},
    {'name': 'double', 'code': 'real(dp)'}
]

with open(argv[1]) as f:
    template = Template(f.read())

out = template.render(types=types)

with open('out.' + argv[1], 'w') as f:
    f.write(out)
