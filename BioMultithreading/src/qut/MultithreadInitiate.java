/**
 * 
 */
package qut;

import java.io.IOException;
import java.util.List;
import java.util.function.Function;

import edu.au.jacobi.pattern.Match;

/**
 * @author aziz
 *
 */
public class MultithreadInitiate implements Runnable {
	// Define variables
	String filename;
	List<Gene> referenceGenes;
	ProgramEntrance getFunctions;
	
	public MultithreadInitiate(String dirString, List<Gene> referenceGenesList) {
		// Initiate variables
		this.filename = dirString;
		this.referenceGenes = referenceGenesList;
		this.getFunctions = new ProgramEntrance();
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
        System.out.println(filename);
        GenbankRecord record = null;
		try {
			record = getFunctions.Parse(filename);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        long startFileTime = System.nanoTime();
        for (Gene referenceGene : referenceGenes)
        {
        	System.out.println(referenceGene.name);
            for (Gene gene : record.genes)
            {
            	Boolean isHomologous = false;
            	isHomologous = getFunctions.Homologous(gene.sequence, referenceGene.sequence);
            	if (isHomologous)
                {
                    NucleotideSequence upStreamRegion = getFunctions.GetUpstreamRegion(record.nucleotides, gene);
                    Match prediction = getFunctions.PredictPromoter(upStreamRegion);
                    if (prediction != null)
                    {
                    	ProgramEntrance.addPrediction(referenceGene, prediction);
                    }
                }
            }
        }
        long eachFileTime = System.nanoTime() - startFileTime;
		System.out.println("Each File: " + (double)eachFileTime / 1000000000.0);
	}
}
