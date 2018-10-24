package qut;

import jaligner.*;
import jaligner.matrix.*;
import edu.au.jacobi.pattern.*;
import java.io.*;
import java.util.*;

public class ProgramEntrance
{
	private Sigma70Definition sigmDef = new Sigma70Definition();
	private BLOSUM62 blosum = new BLOSUM62();
    public static HashMap<String, Sigma70Consensus> consensus = new HashMap<String, Sigma70Consensus>();
    private Series sigma70_pattern = sigmDef.getSeriesAll_Unanchored(0.7);
    private final Matrix BLOSUM_62 = blosum.Load();
    private static byte[] complement = new byte['z'];

    static
    {
        complement['C'] = 'G'; complement['c'] = 'g';
        complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a';
        complement['A'] = 'T'; complement['a'] = 't';
    }

    BioPatterns avoidConflict = new BioPatterns();
    SmithWatermanGotoh avoidConflictGotoh = new SmithWatermanGotoh();
    
    private static List<Gene> ParseReferenceGenes(String referenceFile) throws FileNotFoundException, IOException
    {
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true)
        {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    public boolean Homologous(PeptideSequence A, PeptideSequence B)
    {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    public NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene)
    {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
           upStreamDistance = gene.location-1;

        if (gene.strand == 1)
            return new NucleotideSequence(java.util.Arrays.copyOfRange(dna.bytes, gene.location-upStreamDistance-1, gene.location-1));
        else
        {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i=0; i<upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart-i]];
            return new NucleotideSequence(result);
        }
    }

    public Match PredictPromoter(NucleotideSequence upStreamRegion)
    {
        return avoidConflict.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    public static void ProcessDir(List<String> list, File dir)
    {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    public static List<String> ListGenbankFiles(String dir)
    {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    public GenbankRecord Parse(String file) throws IOException
    {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }
    
    public static synchronized void addPrediction(Gene referenceGene, Match prediction) {
    	consensus.get(referenceGene.name).addMatch(prediction);
    	consensus.get("all").addMatch(prediction);
    }
    
    private static void multithreadWork(String dirString, List<Gene> referenceGenesList) throws InterruptedException {
    	List<String> listOfFiles = ListGenbankFiles(dirString);
        int filesCount = listOfFiles.size();
        
        // Create a thread and runnable for each file
        Thread[] threadsToWork = new Thread[filesCount];
        MultithreadInitiate[] runablies = new MultithreadInitiate[filesCount];
        
        /*Initiate the multi-threading process
        distribute the work load*/
        for (int t = 0; t < filesCount; t++) {
			runablies[t] = new MultithreadInitiate(listOfFiles.get(t), referenceGenesList);
			threadsToWork[t] = new Thread(runablies[t]);
			threadsToWork[t].start();
			System.out.println("thread " + t + " is working...");
		}
        
        // Wait for all thread to finish then print the output
        for (int results = 0; results < filesCount; results++) {
        	threadsToWork[results].join();
        	if(filesCount == results + 1) {
        		for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
        	           System.out.println(entry.getKey() + " " + entry.getValue());
        	}
        }
    }

    public static void run(String referenceFile, String dir) throws FileNotFoundException, IOException, InterruptedException
    {
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        
        // Get the amount of processors and print it out
        int avaliableProcessors = Runtime.getRuntime().availableProcessors();
        System.out.println("Processors avaliable on this machine= " + avaliableProcessors);
        
        // Start multi-threading
        multithreadWork(dir, referenceGenes);
    }

    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException
    {
    	long startSeqTime = System.nanoTime();
		run("./referenceGenes.list", "./Ecoli");
        long sequentialTime = System.nanoTime() - startSeqTime;
		System.out.println("Final time: " + (double)sequentialTime / 1000000000.0);
    }
}
