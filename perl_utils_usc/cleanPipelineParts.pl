#!/usr/bin/perl
die "must be a dir" unless -d $ARGV[0];
chdir($ARGV[0]);
@files = `find`;
chomp @files;

@del = grep /s_.+?sequence\.\d+\./, @files;
@del  = (@del, grep(/\d+\..*?fastq$/, @files));
@del  = grep !/tophat/, @del;
@del  = grep !/unmapped/, @del;
print("$_\n") for @del;
system("rm $_") for @del;