def parse_cds_header(header_str):
	data = header_str.split('|')
	data_count = len(data)
	data_dict = dict(zip(data[0:data_count:2], data[1:data_count:2]))
	return data_dict

def parse_nt_header(header_str):
	data = header_str.split('|')
	data_dict = {}
	data_dict['gi'] = int(data[1])
	data_dict['accession'] = data[3]
	data_dict['description'] = data[4].strip()
	return data_dict

def parse_genome_header(header_str):
	hdata = header_str.split()
	data = parse_cds_header(hdata[0])
	if len(hdata) == 2:
		data['description'] = hdata[1]
	else:
		data['description'] = ''
	return data

def get_cds_id(cds_header_data):
	return '{0}_{1}'.format(cds_header_data['gb'], cds_header_data['cds'])
