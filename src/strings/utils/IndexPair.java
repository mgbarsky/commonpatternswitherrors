package strings.utils;

import java.io.Serializable;

public class IndexPair implements Serializable
{
/**
 *Simple class consists of two integer fields. Can be used in any place where two-dimensional point coordinates or 
 *the start and and end indices of a substring should be saved or passed.
 *The class implements Serializable - in order to be able to save and to read from secondary storage.
 */
	private static final long serialVersionUID = -3230705850675333275L;
	private int _index1;
	private int _index2;		

	public IndexPair(int index1, int index2)	{
		_index1=index1;	
		_index2=index2;
	}

	public int getIndex1()	{	return _index1;	}

	public int getIndex2()	{	return _index2;	}

	public String toString() {	return "("+_index1+","+_index2+")";	}
}


