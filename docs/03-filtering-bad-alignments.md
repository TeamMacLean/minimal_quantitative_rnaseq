# Filtering Badly Aligned Reads

Once we have an alignment, the next step is often to throw out the reads that align badly or not in pairs as we we expect. To do this we need to look at the alignments and assess them one-by-one. We'll need first to have some understanding of the output from our alignment, in this case `aln.sam` a SAM format file.

## SAM Format

Alignments are generally stored in SAM format, a standard for describing how each read aligned one-by-one. Each line carries the results for a single read. Let's examine a single reads alignment. Recall that we can look at one line in a file called `aln.sam` using `tail -n 1 aln.sam` (this gives the bottom line in the file). Running this prints the following

```
NC_011750.1_1004492_1005000_1:0:0_3:0:0_1869f	147	NC_011750.1	1004931	33	70M	=	1004492	-509 TTATATTATTTGGGTTCCTGTGCTGGCGGCTATCTGGAGTATTGGCAGCCTGACAAGCAATGCCTACAAA 2222222222222222222222222222222222222222222222222222222222222222222222	NM:i:3	ms:i:110	AS:i:110	nn:i:0 tp:A:P	cm:i:2	s1:i:59	s2:i:0	de:f:0.0429	rl:i:0
```

On close inspection we can see this mess (which is only a single line) contains things like the read name, the position it maps to on the reference sequence, the read sequence, and lots of other strange things like `70M` and `de:f:0.0429`. The important thing to note is that these weird things are encoded quality information for this alignment, so we can - if we know how to manipulate those codes - select read alignment of the proper quality.

Thankfully the program `samtools` makes this easy for us.

## `samtools`

We can accomplish read filtering with the following command. 

```
samtools view -S -h -q 25 -f 3 aln.sam > aln.filtered.sam

```
Try running that and looking at the output file that is generated. You should have another SAM format file called `aln.filtered.sam` in your working directory. 

Let's take a look at that command in detail

## The `samtools` command and options

Straight away, the command seems to fit the familiar `program name` `options` `files` pattern. It starts with 

```
samtools
```

which is the program name. Then we get the options

```
          view -S -h -q 25 -f 3
```

The first option to `samtools` must be the name of the sub-program to run. There are lots of these as `samtools` is a suite of sub-programs. `view` is the option for working with alignments directly. The second option `-S` tells `samtools view` that we are handing it a SAM format file (soon we will hand it a different type) and `-h` tells it to show the header as well (each SAM file has a header that we sometimes don't want). The next two options are the important ones. `-q 25` will remove reads with a mapping quality (a measure of how well a read is aligned) lower than 25 (a reasonable score) and `-f 3` is a 'flag' a really complex way of encoding alignment attributes (see Further Reading for more details). The important thing is that `3` means `keep reads that are paired and whose pair is mapped too`.   


At the end of the command is the input and output file information

```
                               aln.sam > aln.filtered.sam
```
which means the input file is our `aln.sam` and that the output should be redirected to `aln.filtered.sam`

## Checking the filtering

As an exercise to show that we did filter stuff out lets compare the input `aln.sam` file with the output `aln.filtered.sam` file. Recall that `wc -l` will give us the number of lines in a text file. Run it like this, on both files at once

```
wc -l aln.sam aln.filtered.sam

```

I get this as output

```
  200002 aln.sam
  166905 aln.filtered.sam
  366907 total
```

The number of lines (alignments) in the filtered files is less than that in the unfiltered, so we can casually assume the command worked.

And that's all there is to getting the reads filtered. In real-life you have many options for filtering and you may choose to do it at other points (for instance, lots of RNAseq quantification programs will allow you to filter when you use them), but the process will be similar and take advantage of the same mapping quality and flag metrics you've been introduced to here.

## Are we done?

On the face of it then, it looks like we've come to the end of what we intended to do - we did an alignment, and we've filtered out the poor ones. In practice though, we'll be dealing with many millions of reads, many files of many Gb size. This complicates the housekeeping we have to do, not the procedure we've learned _per se_, so before we jump to the HPC we need to look at that. That's the next chapter.


## Further Reading

### SAM Format

I only really alluded to the SAM format above, but there's a lot to it. This [Wikipedia page](https://en.wikipedia.org/wiki/SAM_(file_format)) gives a lot of detail. 

### Mapping Quality

A metric that describes how well overall the read aligned, it takes into account not just the alignment, but the nubmer of other possible alignments that were rejected. Consider that a read mapping well equally at a number of places in the genome cannot be said to be mapping well at all. Different aligners make arbitrary decisions about how to score such alignments. See this short [summary](https://genome.sph.umich.edu/wiki/Mapping_Quality_Scores) for information on how it can be calculated.
  
### Flags

The flags option is the most powerful way to describe a filter to `samtools view`, it is also really complicated. The number you pass (e.g `-f 3`) is calculated as a sum of lots of options. The way they're are described in the [documentation](https://en.wikipedia.org/wiki/SAM_(file_format)#Bitwise_Flags) is a bit more complex than I want to go into, but there are helpful web-apps that can simplify things - [try this one](https://broadinstitute.github.io/picard/explain-flags.html) 


