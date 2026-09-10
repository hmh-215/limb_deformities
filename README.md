<h1>Whole exome sequencing analyses for detection of limb deformities in human</h1>

<h2> Background </h2>
	<p class="subtitle">Polydactyly and syndactyly are among the most common congenital defects, with the appearance of more than five or the fusion of digits of the hands and feet. Although polydactyly and syndactyly phenotypes are studied extensively, it remains unclear of the underlying genetic factors. Up until now, we were able to take a peak of some pathways involved in the embryonic limp development process, such as <i>Bmp</i>, <i>Wnt</i>, <i>hh</i>, and many others. Remarkably, many of the genetic variants determined as causative for polydactyly and syndactyly mostly affect the <i>hh</i> pathway.</p>
	<p>The objective of these workflows are to perform variant calling for subsequent screening of candidate variants causing limb deformities.</p>

<h2 id="conda-envs">Conda environments</h2>
	<table>
		<tr><th>Environment</th><th>Key tools</th></tr>
		<tr><td><code>preprocessing</code></td><td>Trimmomatic, FastQC, MultiQC</td></tr>
		<tr><td><code>mapping</code></td><td>BWA-MEM, SAMtools, Picard</td></tr>
		<tr><td><code>assembly</code></td><td>SPAdes, QUAST, seqkit, Unicycler</td></tr>
		<tr><td><code>Hcalling</code></td><td>GATK4, Delly, bcftools (+ tabix/bgzip), cnvkit</td></tr>
	</table>
	<p>ANNOVAR (<code>table_annovar.pl</code>, <code>convert2annovar.pl</code>, <code>annotate_variation.pl</code>) is not on Bioconda and is called directly via <code>perl</code> from a manually installed copy under <code>tools/annovar/</code>.</p>
	<p>For unspecified baits as in <code>syndactyly_7_patients.bash</code> pipeline, the script <code>guess_bait.py</code> is used for CNVs calling. Scripts are from <a href="https://github.com/etal/cnvkit">cnvkit repository</a></p>

<h2 id="repo-structure">Pipelines</h2>
	<table>
		<tr><th>Script name</th><th>Purpose</th></tr>
		<tr><td><code>bash/polydactyly/III_2aHNWX.bash</code></td><td>Preprocessing through BQSR, variant calling, ANNOVAR annotation for a Vietnamese polydactyly type II/III patient </td></tr>
		<tr><td><code>bash/syndactyly/s07c.bash</code></td><td>Preprocessing through BQSR, variant calling, ANNOVAR annotation for a Vietnamese syndactyly patient</td></tr>
		<tr><td><code>bash/syndactyly/syndactyly_7_patients.bash</code></td><td>Preprocessing, short-variant / SV / CNV calling and annotation of a 7 syndactyly patients cohort</td></tr>
	</table>

<div class="card">
		<div class="card-title">
			<h3><code>bash/polydactyly/III_2aHNWX.bash</code></h3>
			<span class="tag">Human &middot; WES</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> <code>&lt;sample_id&gt;_1.fastq.gz</code> / <code>_2.fastq.gz</code></span>
			<span><strong>Reference:</strong> hg38</span>
		</div>
		<p>Single-sample pipeline: FastQC/Trimmomatic &rarr; BWA-MEM &rarr; Picard MarkDuplicates &rarr; GATK BaseRecalibrator/ApplyBQSR &rarr; HaplotypeCaller &rarr; SNP/indel splitting and hard-filtering &rarr; ANNOVAR annotation (refGene, cytoBand, ExAC, avsnp150, dbNSFP). Reference indexing and ANNOVAR database downloads are one-time, guarded steps.</p>
		<pre><code>bash ./bash/polydactyly/III_2aHNWX.bash</code></pre>
	</div>

<div class="card">
		<div class="card-title">
			<h3><code>bash/syndactyly/s07c.bash</code></h3>
			<span class="tag">Human &middot; WES</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> <code>&lt;sample_id&gt;_1.fastq.gz</code> / <code>_2.fastq.gz</code></span>
			<span><strong>Reference:</strong> hg38</span>
		</div>
		<p>Single-sample pipeline - with similar structure to <code>bash ./bash/polydactyly/III_2aHNWX.bash</code></p>
		<pre><code>bash ./bash/syndactyly/s07c.bash</code></pre>
	</div>

<div class="card">
		<div class="card-title">
			<h3><code>bash/syndactyly/syndactyly_7_patients.bash</code></h3>
			<span class="tag">Human &middot; WES cohort</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> <code>S&lt;id&gt;.R1/R2.fastq.gz</code> for 7 patients</span>
			<span><strong>Reference:</strong> hg38</span>
		</div>
		<p>Combines preprocessing (BWA-MEM &rarr; sort &rarr; fixmate &rarr; read groups &rarr; dedup, producing separate MQ&ge;10 and MQ&ge;1 BAMs) with three parallel calling arms:</p>
		<ul>
			<li><strong>Short variants:</strong> GATK HaplotypeCaller &rarr; SelectVariants &rarr; hard filtering &rarr; per-cohort merge &rarr; ANNOVAR (refGeneWithVer, cytoBand, gnomAD, avsnp151, dbNSFP, ClinVar).</li>
			<li><strong>Structural variants:</strong> Delly call &rarr; merge &rarr; genotype &rarr; cohort merge &rarr; germline filter &rarr; ANNOVAR (+ DGV).</li>
			<li><strong>CNVs:</strong> CNVkit batch (with optional bait-region inference) &rarr; segmetrics &rarr; call &rarr; per-sample VCF export &rarr; cohort merge &rarr; ANNOVAR (+ DGV).</li>
		</ul>
		<pre><code>bash ./bash/syndactyly/syndactyly_7_patients.bash</code></pre>
	</div>

<h2>Publications</h2>
<ul>
	<li>Nguyen, T. N., & Huynh, M. H. (2024). Comparison of Galaxy and Unix tools for analyzing the exome sequencing data from syndactyly abnormalities. Vietnam Journal of Science and Technology. <a href="https://vjst.vast.vn/jst/article/view/20054">https://doi.org/10.15625/2525-2518/20054</a></li>
</ul>
	
<h2>Notes</h2>
<footer>
		Internal lab pipelines &middot; run on <code>/storage/student9/</code>.
</footer>

</html>
