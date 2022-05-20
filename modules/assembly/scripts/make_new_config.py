import json, os

config = snakemake.params.config
default = snakemake.params.default
outfile = snakemake.output[0]

jdata = {}
cwd = os.getcwd()

for sample, sdata in config.items():

	assembler = config[sample].get('assembler', default)

	jdata[sample] = {
		'r1': config[sample]['r1'],
		'r2': config[sample]['r2'],
		'assembly': os.path.join(cwd, 'assembly', 'completed', assembler, sample + '.fa')
	}

with open(outfile, 'w') as f:
	json.dump(jdata, f, indent=4, sort_keys=True)