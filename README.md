# longread_bcftools_filter
making bcftools filtering easy. implementation of the bcftools which will allow for the faster filtering of the variant calls according to the allelic depth and the tags using simple to overlap approaches as compared to implementing the regular patterns. a bcf filter for filtering the variants and it parses the variants in a flash of a minute. you can filter the variants according to the svtype, filter_type, position, allelicdepth, readdepth and chromosome. A complete human genome parsed in less than 1 minute the link to the data is also given. This is a demo dataset. I am releasing a new faster scala implementation which will be more faster than this and can parse the complete in a flick of a second.

```
bcffiler(bcf_file="/Users/gauravsablok/Desktop/CodeCheck/broad/hs37d5.HG002-SequelII-CCS.bnd-only.sv.vcf", 
                                           svtype = "DUP",readdepth = "21")
1	10448	pbsv.BND.1:10448-10:135524737	C	C]10:135524737]	.	NearReferenceGap	SVTYPE=INV;CIPOS=-9,20;MATEID=pbsv.BND.10:135524737-1:10448	GT:AD:DP	0/1:27,3:30	GT	AD	DP	structural_type	CIPOS	mateID
1	1	1584767.0	pbsv.BND.1:1584767-hs37d5:12063285	C	C[hs37d5:12063285[	.	Decoy	SVTYPE=DUP;CIPOS=0,1;MATEID=pbsv.BND.hs37d5:12...	GT:AD:DP	0/1:16,5:21	[0, 1]	16,5	21	DUP	[CIPOS, 0,1]	[MATEID, pbsv.BND.hs37d5:12063285-1:1584767]
2	1	2368006.0	pbsv.BND.1:2368006-hs37d5:6505668	C	C]hs37d5:6505668]	.	Decoy	SVTYPE=DUP;CIPOS=-35,36;MATEID=pbsv.BND.hs37d5...	GT:AD:DP	0/1:9,2:11	[0, 1]	9,2	11	DUP	[CIPOS, -35,36]	[MATEID, pbsv.BND.hs37d5:6505668-1:2368006]
bcffiler(bcf_file="/Users/gauravsablok/Desktop/CodeCheck/broad/hs37d5.HG002-SequelII-CCS.bnd-only.sv.vcf", svtype = "DUP")
	1	10448	pbsv.BND.1:10448-10:135524737	C	C]10:135524737]	.	NearReferenceGap	SVTYPE=INV;CIPOS=-9,20;MATEID=pbsv.BND.10:135524737-1:10448	GT:AD:DP	0/1:27,3:30	GT	AD	DP	structural_type	CIPOS	mateID
1	1	1584767.0	pbsv.BND.1:1584767-hs37d5:12063285	C	C[hs37d5:12063285[	.	Decoy	SVTYPE=DUP;CIPOS=0,1;MATEID=pbsv.BND.hs37d5:12...	GT:AD:DP	0/1:16,5:21	[0, 1]	16,5	21	DUP	[CIPOS, 0,1]	[MATEID, pbsv.BND.hs37d5:12063285-1:1584767]
2	1	2368006.0	pbsv.BND.1:2368006-hs37d5:6505668	C	C]hs37d5:6505668]	.	Decoy	SVTYPE=DUP;CIPOS=-35,36;MATEID=pbsv.BND.hs37d5...	GT:AD:DP	0/1:9,2:11	[0, 1]	9,2	11	DUP	[CIPOS, -35,36]	[MATEID, pbsv.BND.hs37d5:6505668-1:2368006]
```

Gaurav Sablok \
ORCID: https://orcid.org/0000-0002-4157-9405 \
WOS: https://www.webofscience.com/wos/author/record/C-5940-2014 \
RubyGems Published: https://rubygems.org/profiles/sablokgaurav \
Python Packages Published : https://pypi.org/user/sablokgaurav/
