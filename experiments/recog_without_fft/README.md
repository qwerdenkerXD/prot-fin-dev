# Recognition without Fast Fourier Transformation

The idea of this experiment is to get an idea of the influence of the fft to the recognition algorithm.

As the magic of recognition happens by the combinatorial hashing,
why not combining the raw values of the kidera factors instead of their ff-transformed peaks?

---
## Methods (`methods/*`)
<table>
    <th>method</th><th>steps</th>
    <tr>
        <td><code>get_aa_vectors(file)</code></td>
        <td>
            <ol>
                <li>read kidera factor values for all amino acids from <code>file</code></li>
                <li>extend value table with columns for symbols representing multiple amino acids, by forming the mean of the corresponding amino acids' vectors</li>
            </ol>
        </td>
    </tr>
    <tr>
        <td><code>create constellation(seq_vec)</code></td>
        <td>
            <ol>
                <li>for each value in each vector in <code>seq_vec</code>, index of the vector and the value to the constellation map (an array of these tuples)</li>
            </ol>
        </td>
    </tr>
    <tr>
        <td><code>create_hashes(constellation_map, prot_id)</code></td>
        <td>
            <ol>
                <li>
                    for each index-value-pair in the map create combinatorial hashes (anker points) with the next 100 pairs of the map if these correspond to one in the next 10 amino acids (10 amino acids * 10 kidera factors = max 100 possible)<br>
                    the hashes are generated by multiplying the values with 100 to get 10 bits integers and then combine them into 32-bit int like: <br><code>(index-diff)-(value_of_other_pair)-(value)</code>
                </li>
                <li>save index and protein id for each hash</li>
            </ol>
        </td>
    </tr>
    <tr>
        <td><code>score_prots(hashes, database)</code></td>
        <td>
            <ol>
                <li>for each hash, collect for each protein its offsets to its occurences in the protein sequence</li>
                <li>for each protein, group the occurences of the hashes/ankerpoints by their offsets</li>
                <li>the group with the biggest size is the score for a protein, as it is the best fitting constellation of the hashes</li>
            </ol>
        </td>
    </tr>
    <tr>
        <td><code>create_db(prot_file)</code></td>
        <td>
            <ol>
                <li>create a database for all proteins in the file by joining the results of <code>create_hashes</code></li>
                <li>create a protein-index-map as well to get to information for each protein</li>
            </ol>
        </td>
    </tr>
    <tr>
        <td><code>find_match(fasta_file)</code></td>
        <td>for each protein in the file, find the best scored match</td>
    </tr>
</table>

### Convenience Methods
`hashes_from_seq(seq, aa_vec_map, prot_id)`
 - just the workflow `seq_to_vectors` $\rightarrow$ `create_constellation` $\rightarrow$ `create_hashes`

`seq_to_vectors(seq, aa_vec_map)`
 - transform the sequence into an array of vectors

`iter_fasta(fasta_file)`
 - a generator function to iterate through the fasta file's contents

---
## Results (`results/*`)
|      file      |     content
|----------------|------------------
|initial_test.out|just a test of the basic functionality of the recognition in a set of 1000 proteins<br>input was 6.5% of the original protein sequence of A0A1D8EJF9<br>the first match is the correct one

---
## Discussion/Brainstorming
The initial test implies a fast recognition of similar sequences. But there are tests necessary to check the results with a bigger database or samples from other databases.

Some ideas for future development:
 - add appropriate testing
 - possible to combine more than two points into one hash (maybe `int(value * 10)` $\rightarrow$ 5bit $\rightarrow$ combine with upto 3 others in 20bit)
 - check how to find similarity by specific kidera factors without recreating the database (maybe just allow combinations between the specifics only)

---
## Environment

System: `Ubuntu 20.04.6 LTS`
Shell: `zsh 5.8`

| dependency | version |
|------------|---------|
|   python3  | 3.8.10  |
|    scipy   | 1.10.1  |
|    numpy   | 1.23.0  |
|   pandas   |  2.0.1  |
