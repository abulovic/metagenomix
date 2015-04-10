from __future__ import absolute_import
import json

from metagenomix import supported_aln_types, supported_db_types

_default_fields = {"input_seq_files": list,
				   "aln_type": str,
				   "host_aln_type": str,
				   "microorganism_aln_type": str,
				   "db_type": str,
				   "microorganism_db_type": str,
				   "host_db_type": str,
				   "host_separated": bool,
				   "alignments": list,
				   "host_alignments": list,
				   "microorganism_alignments": list,
				   "supposed_host_tax": int,
				   "output_dir": str}

_mandatory_fields = ("input_seq_files", "host_separated", "output_dir")

_field_missing_str = "Mandatory field '%s' missing."


class ConfigParsingError(Exception):
	pass

def generate_empty_config_file(config_file):
	with open(config_file, "w") as fout:
		config = {key: val() for key, val in _default_fields.iteritems()}
		file_str = json.dumps(config, indent=4)
		fout.write(file_str)

def parse_config(config_file):
	with open(config_file, 'r') as fin:
		try:
			config = json.load(fin)
		except ValueError, e:
			raise ConfigParsingError("Error loading config file,\
				check if file syntax is corrent <%s>" %  config_file)
	for field in _mandatory_fields:
		if not config[field]:
			raise ConfigParsingError(_field_missing_str % field)
	if config["host_separated"]:
		assert(config["host_db_type"] in supported_db_types)
		assert(config["microorganism_db_type"] in supported_db_types)
		for field in ("host_db_type", "microorganism_db_type",
					  "host_alignments", "microorganism_alignments",
					  "host_aln_type", "microorganism_aln_type"):
			if not config[field]:
		   		raise ConfigParsingError("host_separated set to True, but missing\
		   			values in %s field." % field)
	else:
		assert(config["db_type"] in supported_db_types)
		for field in ("db_type", "alignments", "aln_type"):
			if not config[field]:
				raise ConfigParsingError("host_separated set to False, but missing\
					values in %s field." % field)
	return config
