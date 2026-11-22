import sys

def parse_samplesheet(fname):
    records = []
    with open(fname, "r") as f:
        header = f.readline().rstrip().split(',')
        for line in f:
            row = line.rstrip().split(',')
            record = {}
            record[header[0]] = row[0]
            record[header[1]] = row[1]
            record[header[2]] = row[2]
            records.append(record)
    return records

def add_meta_info(records):
    for record in records:
        if record["fastq_1"] and record["fastq_2"]:
            record["single_end"] = '0'
        elif record["fastq_1"] and not record["fastq_2"]:
            record["single_end"] = '1'

def write_sample_sheet(records, fname):
    rows = []
    headers = ["sample", "single_end", "fastq_1", "fastq_2"]
    rows.append(','.join(headers))
    for record in records:
        row = [record[headers[i]] for i in range(0,len(headers))]
        rows.append(','.join(row))
    with open(fname, "w") as f:
        f.write('\n'.join(rows)+'\n')

if __name__ == "__main__":
    sample_sheet = sys.argv[1]
    out_sample_sheet = sys.argv[2]
    data = parse_samplesheet(sample_sheet)
    add_meta_info(data)
    write_sample_sheet(data, out_sample_sheet)
