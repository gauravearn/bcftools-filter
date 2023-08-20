import pandas as pd
def bcffiler(bcf_file=None, svtype=None, \
                               filter=None, \
                                    chromosome=None, \
                                           position=None, \
                                                allelicdepth=None, \
                                                       readdepth=None):
    """summary_line
    a bcffilter for filtering the variants and it parses the variants in a 
    flash of a minute. you can filter the variants according to the svtype,
    filter_type, position, allelicdepth, readdepth and chromosome. 
    
    Keyword arguments:
    argument -- description
    svtype: type of structual type
    filter: type of ["Decoy", "NearReferenceGap","NearContigEnd", \
                                        "InsufficientStrandEvidence", "NotFullySpanned"]
    position: any specific position you are interested
    allelicdepth: range of allelicdepth
    readdepth: the readdepth value                                   
    Return: return_description
    a writable dataframe and also a vcf file
    """
    
    with open(bcf_file, "r") as openvcf:
        with open(bcf_file + "change", "w") as writevcf:
            writevcf.write("CHROM" + "\t" + "POS"  + "\t" + "ID"  "\t" + "REF" + 
                                         "\t" +  "ALT" +  "\t" + "QUAL" + "\t" + "FILTER" + 
                                               "\t" +"INFO"  + "\t" + "FORMAT" + "\t" +"HG002""\n")
            for line in openvcf.readlines():
                if line.startswith("#"):
                    continue
                writevcf.write(line)
    read_vcf = pd.read_csv(bcf_file + "change", sep = "\t")
    read_vcf["GT"] = read_vcf["HG002"].apply(lambda n: n.split(":")[0].split("/"))
    read_vcf["AD"] = read_vcf["HG002"].apply(lambda n: n.split(":")[1])
    read_vcf["DP"] = read_vcf["HG002"].apply(lambda n: n.split(":")[2])  
    read_vcf["structural_type"] = read_vcf["INFO"].apply(lambda n: \
                                           str(n.split(";")[0])). \
                                                          apply(lambda n: n.replace("SVTYPE=", ""))
    read_vcf["CIPOS"] = read_vcf["INFO"].apply(lambda n: n.split(";")[1]). \
                                                                 apply(lambda n: n.split("="))
    read_vcf["mateID"] = read_vcf["INFO"].apply(lambda n: n.split(";")[2]). \
                                                                    apply(lambda n: n.split("="))
    if bcf_file and svtype:
        svtype_read = ["BND", "INV", "DUP", "CNV"]
        svtype_selected = ''.join([i for i in svtype_read if i == svtype])
        selected_read_vcf = read_vcf.where(read_vcf["structural_type"] == svtype_selected).dropna()
        return selected_read_vcf
    if bcf_file and filter:
        filter_read = ["Decoy", "NearReferenceGap","NearContigEnd", \
                                        "InsufficientStrandEvidence", "NotFullySpanned"]
        filter_store = ''.join([i for i in filter_read if i == filter])
        selected_filter_vcf = read_vcf.where(read_vcf["FILTER"] == filter_store).dropna()
        return selected_filter_vcf
    if bcf_file and position:
       store_position = int(position)
       selected_position = read_vcf.iloc[::].where(read_vcf["POS"] == store_position).dropna()
       return selected_position
    if bcf_file and allelicdepth:
        store_allelic_depth = allelicdepth
        selected_allelicdepth = read_vcf.where(read_vcf["AD"].apply(lambda n: \
                                                            n[0] == store_allelic_depth)).dropna()
        return selected_allelicdepth
    if bcf_file and readdepth:
        store_readdepth = readdepth
        return read_vcf.iloc[::].where(read_vcf["DP"] == store_readdepth).dropna()
    if bcf_file and svtype and readdepth:
        svtype_selected = ''.join([i for i in svtype_read if i == svtype])
        store_read_depth = readdepth
        return read_vcf.iloc[::].where((read_vcf["structural_type"] == svtype_selected) \
                                               .dropna() & (read_vcf["DP"] == store_read_depth)).dropna()
    if bcf_file and svtype and allelicdepth:
        store_allelic_depth = allelicdepth
        svtype_selected = [i for i in svtype_read if i == svtype]
        return read_vcf.where(read_vcf["structural_type"] == svtype_selected).dropna() \
                                                    & (read_vcf["AD"].apply(lambda n: n[0] \
                                                                        == store_allelic_depth)).dropna()
