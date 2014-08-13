commonpatternswitherrors
========================
<h1>Approximate common patterns</h1>
<p>This is an implementation of the algorithm 
for finding in two input texts all pairs of approximate common patterns, 
which satisfy the following user-specified criteria: <br>
<ul>
<li>the length of each pattern is at least <em>s</em> characters. </li>
<li>each pattern occurs in both input texts, but with up to <em>k</em> errors 
(mismatches, insertions or deletions)</li>
</ul>
</p>

<p>The algorithm - All Paths Below the Threshold (APBT) - is described in <br>
<em>M. Barsky, U. Stege, A. Thomo, C. Upton.</em><br> 
<strong>A graph approach to the threshold all-against-all substring matching problem.</strong><br>
ACM Journal of Experimental Algorithmics 12: 2008.
</p>

<p>This is an in-memory algorithm with the running time <em>O(NMk^3)</em>, 
where <em>N</em> and <em>M</em> are the lengths of two input texts. 
Thus its scalability is limited. However, this is the fastest algorithm for this problem, 
as the previous algorithm was running in <em>O(N^2M^2)</em> time.</p> 

<p>The output consists of pairs of approximate common patterns with their starting positions in each input text.</p>

<p>Program works with texts over any alphabet. It is faster for larger alphabets: for example it is
an order of magnitude faster for protein sequences, than for DNA sequences.</p>

<h2>To compile:</h2>
Create <em>bin</em> folder. From <em>src</em> folder:
<pre><code>
set path=<em>path to JDK</em>/bin
javac -cp . strings/algorithms/APBT.java -d ../bin
</code></pre>
The program was developed with an old version of Java, without generics. 
Ignore compilation warnings.

<h2>Program parameters</h2>
Specify command-line arguments in the following order:
<ol>
<li>Name of the first input file</li>
<li>Name of the second input file</li>
<li>Minimum pattern length <em>s</em></li>
<li>Maximum number of errors <em>k</em></li>
<li>Compress output into maximal patterns? 1-yes, 0-no</li>
<li>Print to stdout? 1-yes, 0-no</li>
</ol>

<h2>To run:</h2>
<pre><code>
java -Xmx512M -Xms512m strings.algorithms.APBT ../sample_inputs/humanprotein.txt ../sample_inputs/mouseprotein.txt 51 3 0 0
</code></pre>
or use the jar-packaged app with the same arguments:
<pre><code>java -Xmx512M -Xms512m -jar APBT.jar arguments</code></pre>
<h2>Sample usage:</h2>
There are two sample DNA files and two sample protein files in a compressed folder 'sample_inputs'. 
Sample runs and sample outputs are recorded in file 'RUN_SAMPLES.txt'
