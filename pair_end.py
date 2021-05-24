import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--datapath', type=str, default='test/', help='Data path')
parser.add_argument('--workpath', type=str, default='work/', help='Work path')
parser.add_argument('--genome', type=str, default='hg19', help='Reference genome')
parser.add_argument('--threads', type=int, default=10, help='Number of threads')
parser.add_argument('--server', type=str, default='../server.txt', help='File for saving the server information')

args = parser.parse_args()
print args

genome = args.genome
thread = args.threads
datapath = args.datapath
workpath = args.workpath
server = {}
with open(args.server, 'r') as tmp:
    for line in tmp:
        ele = line.split()
        server[ele[0]] = ele[1]

os.environ['PATH'] = os.environ['PATH']+':'+os.path.abspath('bin')
os.environ['BOWTIE2_INDEXES'] = '/data/BOWTIE2_INDEXES/'

htmlmain = open(workpath+'index.html', 'w')

with open('../styles.css', 'r') as temp:
    lines = temp.readlines()
htmlmain.write('''
<head>
<style type="text/css">
%s
</style>
'''%''.join(lines))

dataset = []
for p in os.listdir(datapath):
    if os.path.isfile(datapath+p):
        continue
    f = [i for i in os.listdir(datapath+p) if i.endswith('.fastq.gz')]
    f.sort()
    print p, f
    if len(f) == 2:
        dataset.append((p, datapath+p+'/'+f[0], datapath+p+'/'+f[1]))
dataset.sort()

htmlmain.write('''
<h1>Report for Single Cell RNA-seq</h1>
<h2>Data Files</h2>
<table>
<tr><th>Name</th><th>Read 1</th><th>Read 2</th></tr>
''')
for d, f1, f2 in dataset:
    htmlmain.write('<tr><th>%s</th><td>%s</td><td>%s</td></tr>'%(d,f1,f2))

htmlmain.write('''
</table>        
        
<h2>Data Quality</h2>
<table>
<tr><th>Name</th><th>Fastx output</th><th>Fastqc output</th></tr>
 ''')
for d, f1, f2 in dataset:
    print '> Check', d, 'to', genome
    out = '%s.%s.fastx.txt'%(d, genome)
    if not os.path.exists(workpath+out):
        os.system('zcat '+f1+' | fastx_quality_stats > '+workpath+out)
        os.system('zcat '+f2+' | fastx_quality_stats >>'+workpath+out)
    f1name = (f1.split('/')[-1]).split('.')[0]+'_fastqc'
    if not os.path.exists(workpath+f1name):
        os.system('fastqc '+f1+' -o '+workpath)
    f2name = (f2.split('/')[-1]).split('.')[0]+'_fastqc'
    if not os.path.exists(workpath+f2name):
        os.system('fastqc '+f2+' -o '+workpath)
    htmlmain.write('<tr><th>%s</th><td><a href="%s">%s</a></td><td><a href="%s/fastqc_report.html">Read 1</a> <a href="%s/fastqc_report.html">Read 2</a></td></tr>'%(d, out, out, f1name, f2name))
htmlmain.write('''
</table>        

<h2>Map reads by Tophat</h2>
<table>
<tr><th>Name</th><th>Command</th><th>Log</th></tr>
''')
bams = []
for d, f1, f2 in dataset:
    print '> Tophat', d, 'to', genome
    out = '%s.%s.tophat.bam'%(d, genome)
    bams.append(workpath+out)
    cmd = 'tophat2 -G refseq_%s.gtf --transcriptome-index=transidx -g 1 -p %s %s '%(genome, thread, genome) + '%s %s'%(f1,f2)
    if not os.path.exists(workpath+out):
        os.system(cmd)
        os.system('mv tophat_out/accepted_hits.bam '+workpath+out)
        os.system('mv tophat_out/junctions.bed '+workpath+out+'.junctions.bed')
        os.system('rm -rf tophat_out')
    if not os.path.exists(workpath+out+'.log'):
        os.system('samtools flagstat '+workpath+out+'>'+workpath+out+'.log')
    htmlmain.write('<tr><th>%s</th><td>%s</td><td><a href="%s.log">Open</a></td></tr>'%(d, cmd, out))
    if not os.path.exists(workpath+out+'.bai'):
        os.system('samtools index '+workpath+out)
htmlmain.write('''
</table>

<h2>Genome browser tracks</h2>
<p>Copy the following tracks:</p>
<table><tr><td>''')
for d, f1, f2 in dataset:
    print '> BAM to BigWig', d
    out = '%s.%s.tophat.bw'%(d, genome)
    htmlmain.write('track type=bigWig db=%s name="%s" visibility=2 bigDataUrl="%s/%s/single-cell-rna-seq/%s"<br/>'%(genome, d, server['address'], server['home'], out))
    if not os.path.exists(workpath+out):
        os.system("bedtools genomecov -bg -split -ibam work/%s.%s.tophat.bam > tmp.bdg"%(d, genome))
        os.system("bedGraphToBigWig tmp.bdg %s.chrom.sizes "%genome+workpath+out)
        os.system("rm tmp.bdg")

htmlmain.write('''</td></tr></table>
<p>into <a href="http://genome.ucsc.edu/cgi-bin/hgCustom?db=%s" target="_black">the UCSC genome browser</a>, or use <a href="%s/internal/browser/cgi-bin/hgCustom?db=%s" target="_black">the local genome browser</a> (user name: xuelab; password: neibu)</p>

<h2>Read counts by transcripts<h2>
'''%(genome, server['address'], genome))

cmd1 = 'featureCounts -t exon -g gene_id -O -M -a refseq_%s.gtf -o %sfeature_counts %s'%(genome, workpath, ' '.join(bams))
cmd2 = 'sed 1d %sfeature_counts | cut -f1,7- > %sfeature_counts_brief.txt'%(workpath,workpath)

if not os.path.exists(workpath+'feature_counts_brief.txt'):
    os.system(cmd1)
    os.system(cmd2)

refseq2gene = {}
if os.path.exists('refseq2gene_hg19.txt'):
    infile = open('refseq2gene_hg19.txt', 'r')
    for line in infile:
       a,b = line.split('\t') 
       refseq2gene[a] = b.strip()
    infile.close()
    print 'Get', len(refseq2gene), 'refseq ids mapped to gene names'

vals = {}
for i in xrange(len(dataset)):
    d, f1, f2 = dataset[i]
    val = {}
    temp = open(workpath+'feature_counts_brief.txt', 'r')
    temp.readline() ## skip header
    for line in temp:
        ele = line.split('\t')
        v = int(ele[i+1])
        if v >0:
            val[refseq2gene.get(ele[0], 'Unknown')+'\t'+ele[0]] = v
    temp.close()
    vals[d] = val

output = open(workpath+'feature_counts_gene.txt', 'w')
genes = val.keys()
genes.sort()
output.write('Refseq\tGene%s\n'%'\t'.join([d for d,f1,f2 in dataset]))
for gene in genes:
    output.write(gene)
    for d,f1,f2 in dataset:
        val = vals[d]
        if gene in val:
            output.write('\t%s'%val[gene])
        else:
            output.write('\t0')
    output.write('\n')
output.close()

htmlmain.write('''
<table>
<tr><td>%s</td></tr>
<tr><td>%s</td></tr>
<tr><td>Output file is <a href="feature_counts_gene.txt">feature_counts_gene.txt</a></td></tr>
</table>
'''%(cmd1, cmd2))

htmlmain.write('''
<h2>Overlap of expressed genes</h2>
<table>
''')
for d1, f1, f2 in dataset:
    htmlmain.write('<tr><th>'+d1+'</th>')
    v1 = vals[d1]
    htmlmain.write('<td>%d expressed</td>'%len(v1))
    for d2, f1, f2 in dataset:
        v2 = vals[d2]
        g1 = set([g for g in v1 if v1[g]>0])
        g2 = set([g for g in v2 if v2[g]>0])
        htmlmain.write('<td>%.f%% in %s</td>'%(100*len(g1&g2)/float(len(g1)), d2))
    htmlmain.write('</tr>')

htmlmain.write('''
</table>

<h1>End of Report</h1>
''')
htmlmain.close()
