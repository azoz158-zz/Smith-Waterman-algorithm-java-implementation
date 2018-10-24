package qut;

import jaligner.*;
import jaligner.matrix.*;
import java.io.*;


public class PeptideSequence 
{
    public byte[] bytes;
    static BLOSUM62 blosum = new BLOSUM62();

    public final static Matrix BLOSUM_62 = blosum.Load();

    public PeptideSequence()
    {
    }
    
    public PeptideSequence(String string)
    {
        bytes = string.getBytes();
    }
    
    public static double Similarity(PeptideSequence A, PeptideSequence B)
    {  
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore();       
    }
    
    @Override
    public String toString()
    {
        return new String(bytes);
    }
}
