# Running `minimap2`

Running `minimap2` takes only one step. Assuming we've already `cd`'d into the directory with the reads and reference we can use this command

```
  minimap2 -ax sr ecoli_genome.fa ecoli_left_R1.fq ecoli_right_R2.fq > aln.sam
```

Try running that and see what happens... You should get an output file in the working directory called `aln.sam`. On my machine this takes just a few seconds to run.

Let's look at the command in detail.

## The `minimap2` command and options

First we get this

```
  minimap2
```

which is the name of the actual program we intend to run, so it isn't surprising that it comes first. The rest of the command are options (sometimes called arguments) telling the program how to behave and what it needs to know. Next up is this

```
           -ax sr
```

which gives option `a` meaning print out SAM format data. And option `x` meaningwe wish to use a preset parameter set. The preset we wish to use comes after `x` and is `sr`, which stands for `short reads` and tells `minimap2` to use settings for short reads against a long genome. Next is this

```
                   ecoli_genome.fa ecoli_left_R1.fq ecoli_right_R2.fq
```
which are the input files in the 'reference' 'left read' 'right read' order. Finally, we have

```
                                                                       > aln.sam
```

which is the `>` output redirect operator and the name of an output file to write to. This bit specifies where the output goes.

So the structure of the `minimap2` command (like many other commands) is simply `program_name options input output`. 

And this one command is all we need for a basic alignment with `minimap2`. We can now move on to the next step in the pipeline.


## Further Reading

### The `>` operator

The `>` symbol is actually not part of the `minimap2` command at all, it is a general shortcut that means something like 'catch the output from the process on the left and put it in the file on the right`. Think of the `>` as being a physical funnel catching the datastream! Because it's a general operator and not an option in a program, we can almost always use '>' to make output files. You'll see it pop up quite often

### `minimap2` further instructions and github

The commands given here for `minimap2` are just a small selection of what are available. You can see the user guide at [GitHub](https://github.com/lh3/minimap2)
