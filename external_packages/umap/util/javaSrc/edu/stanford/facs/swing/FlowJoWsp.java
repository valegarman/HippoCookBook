/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 */
package edu.stanford.facs.swing;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Properties;
import java.util.TreeMap;
import java.util.List;

import org.apache.xerces.dom.*;
import org.w3c.dom.Element;

public  class FlowJoWsp {
	
	public final Properties props, propsGui;

	public static final String 
			SAMPLE="Sample",
			FOCI="gating:foci",
			EDGE="gating:edge", 
			VERTEX="gating:vertex", 
			COORD="gating:coordinate",
			TYPE_ROOT="root", 
			TYPE_FOLDER="folder",  
			TYPE_SAMPLE="sample", 
			TYPE_GATE="gate",
			ROOT_ID=TYPE_ROOT+":0",
			ROOT_CHILDREN=ROOT_ID+".children";
	        
	public boolean IsGateId(final String id) {
		return id.startsWith(TYPE_GATE +':');
	}

	public boolean IsSampleId(final String id) {
		return id.startsWith(TYPE_SAMPLE  +':');
	}

	public boolean IsFolderId(final String id) {
		return id.startsWith(TYPE_FOLDER  +':');
	}

	public boolean IsRootId(final String id) {
		return id.startsWith(TYPE_ROOT  +':');
	}

	public static final String []
			ELLIPSE_FOCI= {FOCI, VERTEX, COORD},
			ELLIPSE_EDGE= {EDGE, VERTEX, COORD};
	public final DeferredDocumentImpl doc;
	public final String version, flowJoVersion;
	public DeepNodeListImpl sampleNodes;
	public int nSamples;
	public final Map<String, ElementImpl> idMap=new TreeMap<String, ElementImpl>();
	
	public void removeChild(final String pid, final String child) {
		final String prop=pid+".children";
		final String str=this.props.getProperty(prop);
		if (!Basics.isEmpty(str))
			this.props.setProperty(prop, str.replace(child+" ", ""));
	}
	
	public void addChild(final String pid, final String childId) {
		final String prop=pid+".children";
		String p=props.getProperty(prop);
		if (p==null)
			p="";
		props.setProperty(prop, p+childId +" ");
		//System.out.println("pid="+pid+", childId="+childId+", prop="+prop+", getProperty(prop)="+props.getProperty(prop));
	}
	
	public void setParent(final Object childIdOrSample, final String pid) {
		final String childId;
		if (childIdOrSample instanceof Integer) {
			childId=(String)this.getSampleIdByNum(((Double)childIdOrSample).intValue());
		}else if (childIdOrSample instanceof Integer) {
			childId=(String)this.getSampleIdByNum((Integer)childIdOrSample);
		} else {
			childId=(String)childIdOrSample;
		}
		final String pProp=childId+".parent";
		final String priorParent=this.props.getProperty(pProp);
		if (!Basics.isEmpty(priorParent))
			removeChild(priorParent, childId);
		addChild(pid, childId);
		this.props.setProperty(pProp, pid);
		rootIdsInSync=false;
	}
	
	public Object getSampleIdByNum(final int num) {
		return TYPE_SAMPLE+":"+getSampleNode(num).getAttribute("sampleID");
	}
	
	public static Object []parseIdTypeNumb(String id) {
		final int idx=id.indexOf(':');
		if (idx<0) {
			return new String[] {id, "?", id};
		}
		final String type=id.substring(0, idx);
		id=id.substring(idx+1);
		final Object [] out;
		if (id.startsWith("ID"))
			out= new Object[]{id, type, id.substring(2)};
		else
			out=new Object[]{id, type, id};
		return out;
	}
	
	public static Object parseId(final String id) {
		final int idx=id.indexOf(':');
		if (idx<0) {
			return id;
		}
		return id.substring(idx+1);
	}
	
	public void refreshSampleNodes() {
		sampleNodes=(DeepNodeListImpl) doc.getElementsByTagName(SAMPLE);
		nSamples=sampleNodes.getLength();
	}
	
	// num is 1 based
	public ElementImpl getSampleNode(final int num) {
		return GetFirstNode((ElementImpl) sampleNodes.item(num-1), "SampleNode");
	}
	
	// return value is 1 based
	public int getSampleNumById(final String sampleId) {
		final Object id=parseId(sampleId);
		for (int i=1;i<=this.nSamples;i++) 
			if (getSampleNode(i).getAttribute("sampleID").equals(id))
				return i;
		return 0;
	}

	public Object get(final Object nodeOrNum, final String []names, final boolean all) {
		final ElementImpl in;
		if (nodeOrNum instanceof Integer) 
			in=(ElementImpl) sampleNodes.item(((Integer)nodeOrNum).intValue()-1);
		else if (nodeOrNum instanceof Double) 
			in=(ElementImpl) sampleNodes.item(((Double)nodeOrNum).intValue()-1);
		else
			in=(ElementImpl)nodeOrNum;
		if (all)
			return GetAllWild(in, names);
		else
			return GetFirst(in, names, 1);
	}
	
	public FlowJoWsp(final DeferredDocumentImpl doc, 
			final Properties props, 
			final Properties propsGui) {
		this.doc=doc;
		this.props=props;
		this.propsGui=propsGui;
		final ElementImpl wsp=(ElementImpl) doc.item(0);
		version=wsp.getAttribute("version");
		flowJoVersion=wsp.getAttribute("flowJoVersion");
		this.refreshSampleNodes();
	}
	
	public static ElementImpl GetFirstNode(final ElementImpl node, 
			final String name) {
		final int N=node.getLength();
		for (int i=0;i<N;i++) {
			if (name.equals(node.item(i).getNodeName())) {
				return (ElementImpl)node.item(i);
			}
		}
		return null;
	}
	
	public static ElementImpl GetFirst(final ElementImpl node, 
			final String []names, final int nameIdx) {
		final String name=names[nameIdx-1];
		final int N=node.getLength();
		for (int i=0;i<N;i++) {
			if (name.equals(node.item(i).getNodeName())) {
				if (nameIdx==names.length)
					return (ElementImpl)node.item(i);
				else
					return GetFirst((ElementImpl) node.item(i), names, nameIdx+1);
			}
		}
		return null;
	}
	
	public static Collection<ElementImpl> GetList(
			final ElementImpl node, 
			final String []names) {
		final Collection<ElementImpl> all=new ArrayList<>();
		GetAll(all, node, names, 1);
		return all;
	}

	public static ElementImpl[]GetAll(
			final ElementImpl node, 
			final String []names) {
		return GetList(node, names).toArray(emptyEle);
	}
	
	public static Collection<ElementImpl> GetAll(
			final Collection<ElementImpl> all,
			final ElementImpl node, 
			final String []names, final int nameIdx) {
		final String name=names[nameIdx-1];
		final int N=node.getLength();
		for (int i=0;i<N;i++) {
			if (name.equals(node.item(i).getNodeName())) {
				if (nameIdx==names.length)
					all.add((ElementImpl)node.item(i));
				else
					GetAll(all, (ElementImpl) node.item(i), names, nameIdx+1);
			}
		}
		return all;
	}

	public static ElementImpl GetFirstWild(
			final ElementImpl node, 
			final String []names, final int nameIdx) {
		String name=names[nameIdx-1];
		final int len=name.length();
		final int N=node.getLength();
		if (len>1 && name.endsWith("*")){
			name=name.substring(0, len-1);
			for (int i=0;i<N;i++) {
				if (node.item(i).getNodeName().startsWith(name)) {
					if (nameIdx==names.length)
						return (ElementImpl)node.item(i);
					else
						return GetFirstWild((ElementImpl) node.item(i), names, nameIdx+1);
				}
			} 
		}else {
			for (int i=0;i<N;i++) {
				if (name.equals(node.item(i).getNodeName())) {
					if (nameIdx==names.length)
						return (ElementImpl)node.item(i);
					else
						return GetFirstWild((ElementImpl) node.item(i), names, nameIdx+1);
				}
			}

		}
		return null;
	}
	private static final ElementImpl []emptyEle=new ElementImpl[0];

	private static final String [] subPop=new String[] {"Subpopulations"};
	
	public Element getSubPopNode(
			final ElementImpl pop) {
		ElementImpl[]s=GetAll(pop, subPop);
		if (s!= null && s.length>0) return s[0];
		final Element node=doc.createElement(subPop[0]);
		pop.appendChild(node);
		return node;
	}
	
	
	private static final String [] graph=new String[] {"Graph"};
	private static final String [] axis=new String[] {"Graph", "Axis"};
	public void addParentGraph(
			final ElementImpl parentPopulation, 
			final ElementImpl childPopulation, 
			final String []dims) {
		
		final ElementImpl []parentGraphs=GetAll(parentPopulation, graph);
		if (parentGraphs != null && parentGraphs.length>0) {
			// duplicate 1st graph ... not sure if FlowJo adds others
			final ElementImpl graph=parentGraphs[0];
			final ElementImpl []axes=GetAll(parentPopulation, axis);
			if (axes==null || axes.length !=2) {
				for (int i=0;i<2;i++) {
					final Element ax=doc.createElement("Axis");
					graph.appendChild(ax);
					ax.setAttribute("auto", "x");
					ax.setAttribute("label", "");
					ax.setAttribute("name", dims[i]);
				}
			} else {
				for (int i=0;i<2;i++) {
					axes[i].setAttribute("name", dims[i]);
				}
			}
			childPopulation.appendChild(parentGraphs[0].cloneNode(true));
		}
	}
	
	public static void setIds(final Element node, 
			final String pid, 
			final String id, 
			final boolean isParentSample) {
		if (!isParentSample) {
			node.setAttribute("gating:parent_id", pid);			
		}
		node.setAttribute("gating:id", id);
	}
	
	
	public ElementImpl setDim(
			final ElementImpl gatingNode, 
			final String dimName) {
		if (gatingNode==null) {
			System.out.print("gatingNode is null");
			return null;
		}
		final ElementImpl dim=(ElementImpl)
				this.doc.createElement("gating:dimension");
		gatingNode.appendChild(dim);
		final ElementImpl fcsDim=(ElementImpl)
				this.doc.createElement("data-type:fcs-dimension");
		fcsDim.setAttribute("data-type:name",  dimName);
		dim.appendChild(fcsDim);
		return dim;
	}

	public static void debugMatrix(double [][]matrix) {
		System.out.println("rows x cols ... " +matrix.length +"x" + matrix[0].length);
	}
	
	public Object[]createPolygonSubGate(
			final ElementImpl parentPopulation, 
			final String parentId,
			final String newGateId, // must be supplied by MATLAB factory
			final double [][]roiPositionUnscaled, 
			final String name,
			final int count, 
			final String []dims){
		final Object []output=new Object[3];
		output[0]=TYPE_GATE+":"+newGateId;
		final ElementImpl newPopulation=(ElementImpl)doc.createElement("Population");
		output[1]=newPopulation;
		newPopulation.setAttribute("name", name);
		newPopulation.setAttribute("count", Integer.toString(count));
		newPopulation.setAttribute("annotation", "");
		newPopulation.setAttribute("owningGroup", "");
		newPopulation.setAttribute("expanded", "0");
		newPopulation.setAttribute("sortPriority", "10");
		getSubPopNode(parentPopulation).appendChild(newPopulation);
		addParentGraph(parentPopulation, newPopulation, dims);
		final ElementImpl newGateNode = (ElementImpl)doc.createElement("Gate");
		output[2]=newGateNode;
		final boolean isSample=parentPopulation.getNodeName().equals("SampleNode");
		final String pid=(String)parseId(parentId);
		this.idMap.put((String)output[0], newPopulation);
		setIds(newGateNode, pid, newGateId, isSample);
		newPopulation.appendChild(newGateNode);
		
		final ElementImpl gatingNode=(ElementImpl)doc.createElement("gating:PolygonGate");
		gatingNode.setAttribute("isTinted", "0");
        gatingNode.setAttribute("eventsInside", "1");
        gatingNode.setAttribute("lineWeight", "Normal");
        gatingNode.setAttribute("tint", "#000000");
        gatingNode.setAttribute("userDefined", "1");
        setIds(gatingNode, pid, newGateId, isSample);
        setDim(gatingNode, dims[0]);
        setDim(gatingNode, dims[1]);
        newGateNode.appendChild(gatingNode);
		gatingNode.setAttribute("quadId", "-1");
        gatingNode.setAttribute("gateResolution", "256");
        // setVertices(gatingNode, roiPosition);
        for (int r=0;r<roiPositionUnscaled.length;r++) {
        	final ElementImpl v=(ElementImpl)doc.createElement("gating:vertex");
    		gatingNode.appendChild(v);
        	for (int c=0;c<roiPositionUnscaled[0].length;c++) {
        		final ElementImpl coord=(ElementImpl)doc.createElement("gating:coordinate");
        		coord.setAttribute("data-type:value", 
        				Double.toString(roiPositionUnscaled[r][c]));
        		v.appendChild(coord);
        	}
        }
		return output;
	}

	public Object[]createPolygonSubGate(final int sampleNum, 
			final String parentId, final String newGateId,  
			final double [][]roiPositionUnscaled, final String name,
			final int count, final String []dims){
		return createPolygonSubGate(getSampleNode(sampleNum),
				parentId, newGateId, roiPositionUnscaled, name, count, dims);
	}
	
	public static List<ElementImpl> GetWildList(
			final ElementImpl node, 
			final String []names) {
		final List<ElementImpl> all=new ArrayList<>();
		GetAllWild(all, node, names, 1);
		return all;
	}
	
	
	public static ElementImpl []GetAllWild(
			final Object node, 
			final String []names) {
		if (!(node instanceof ElementImpl)) return null;
		return GetWildList((ElementImpl)node, names).toArray(emptyEle);
	}	

	static Collection<ElementImpl> GetAllWild(
			final Collection<ElementImpl> all,
			final ElementImpl node, 
			final String []names, final int nameIdx) {
		String name=names[nameIdx-1];
		final int len=name.length();
		final int N=node.getLength();
		if (len>1 && name.endsWith("*")){
			name=name.substring(0, len-1);
			for (int i=0;i<N;i++) {
				if (node.item(i).getNodeName().startsWith(name)) {
					if (nameIdx==names.length)
						all.add((ElementImpl)node.item(i));
					else
						GetAllWild(all, (ElementImpl) node.item(i), names, nameIdx+1);
				}
			} 
		}else {
			for (int i=0;i<N;i++) {
				if (name.equals(node.item(i).getNodeName())) {
					if (nameIdx==names.length)
						all.add((ElementImpl)node.item(i));
					else
						GetAllWild(all, (ElementImpl) node.item(i), names, nameIdx+1);
				}
			}
		}
		return all;
	}

	public Boolean hasNotProps=null;
	public static final String [] SUBPOPS={"Subpopulations", "Population"};
	public static final String [] SUBNOTS={"Subpopulations", "NotNode"};
	
	
	public static boolean IsNotGate(ElementImpl population) {
		return population!=null && population.getNodeName().equals("NotNode");
	}
	
	public Collection<ElementImpl>getSubPopulationList(final ElementImpl node){
		final Collection<ElementImpl>out=GetList(node, SUBPOPS);
		if (hasNotProps==null) {
			hasNotProps=doc.getElementsByTagName("NotNode").getLength()>0;
		}
		if (this.hasNotProps) {
			out.addAll(GetList(node, SUBNOTS));
		}
		return out;
	}
	
	public ElementImpl[]getSubPopulations(final ElementImpl node){
		return getSubPopulationList(node).toArray(emptyEle);
	}
	
	
	public String ensureUniquePopulationName(final ElementImpl population, final String name) {
		final Collection<ElementImpl> l=getSubPopulationList(population);
		for (final ElementImpl pop:l) {
			if (pop.getAttribute("name").equals(name)) {
				int suffix=2;
				while(true) {
					boolean found=false;
					final String tryName=name+"#"+suffix;
					for (final ElementImpl pop2:l) {
						if (pop2.getAttribute("name").equals(tryName)) {
							found=true;
							break;
						}
					}
					if (!found) {
						return tryName;
					}
					suffix++;
				}
			}
		}
		return name;
	}
	
	public ElementImpl findPopulation(
			final ElementImpl population, 
			final String []names, 
			final int level/*1 based*/) {
		final Collection<ElementImpl> l=getSubPopulationList(population);
		String name=names[level-1];
		if (name.endsWith("*")) {
			name=name.substring(0,name.length()-1);
			for (final ElementImpl pop:l) {
				if (pop.getAttribute("name").startsWith(name)) {
					if (level==names.length)
						return pop;
					else
						return findPopulation(pop, names, level+1);
				}
			} 
		}else {
			for (final ElementImpl pop:l) {
				if (pop.getAttribute("name").equals(name)) {
					if (level==names.length)
						return pop;
					else
						return findPopulation(pop, names, level+1);
				}
			}
		}	
		return null;
	}
	
	public ElementImpl getPopulationById(final String id) {
		final Object []s=parseIdTypeNumb(id);
		final DeepNodeListImpl gates= (DeepNodeListImpl)doc.getElementsByTagName("Gate");
		final int N=gates.getLength();
		for (int i=0;i<N;i++) {
			final ElementImpl gate=(ElementImpl)gates.item(i);
			final String id2=gate.getAttribute("gating:id");
			final ElementImpl pop=(ElementImpl)gate.getParentNode();
			if (s[0].equals(id2))
				return pop;
			idMap.put(s[1]+":"+id2, pop);
		}
		return null;
	}
	
	public ElementImpl getNodeById(final String id) {
		ElementImpl node=idMap.get(id);
		if (node==null) {
			if (IsSampleId(id))
				node=getSampleNode(getSampleNumById(id));
			else if (IsGateId(id))
				node=getPopulationById(id);
			if (node != null)
				idMap.put(id, node);
		}
		return node;
	}
	
	public Object []getNodeAndTypeById(final String id){
		final ElementImpl node=getNodeById(id);
		return new Object[] {node, getGateType(node)};
	}
	
	public Object getId(final ElementImpl node) {
		if (node.getNodeName().equals("SampleNode")) 
			return TYPE_SAMPLE+":"+node.getAttribute("sampleID");
		else
			return getGateId(node);
	}
	
	public Object getParentId(final String id) {
		String pid=props.getProperty(id+".parent");
		if (Basics.isEmpty(pid)) {
			if (!IsSampleId(id) && !IsFolderId(id)) {
				ElementImpl node=getNodeById(id);
				if (node == null)return "";
				pid=(String)getId((ElementImpl)
					node.getParentNode());
				while(Basics.isEmpty(pid)) {
					node=(ElementImpl)node.getParentNode();
					if (node==null)
						break;
					pid=(String)getId((ElementImpl)node);
				}
				props.setProperty(id+".parent", pid);
			}
		}
		return pid;
	}
	
	boolean addGateOrSampleId(final ElementImpl node, final Collection<String>c, final boolean sampleToo) {
		final boolean found="SampleNode".equals(node.getNodeName());
		if (found) {
			if (sampleToo)
				c.add(TYPE_SAMPLE+":"+node.getAttribute("sampleID"));
		}else {
			final String id=(String)getGateId(node);
			if (!Basics.isEmpty(id)) {
				c.add(id);
			}
		}
		return found;
	}
	
	
	public ArrayList<String>getParentIds(ElementImpl node, final boolean sampleToo, final boolean folderToo){
		final ArrayList<String>c=new ArrayList<>();
		boolean found=addGateOrSampleId(node, c, sampleToo);
		if (c.size()>0) {
			String pid=c.get(0);
			if (this.props.containsKey(pid+".parent")) {
				pid=(String)getParentId(pid);
				while (!Basics.isEmpty(pid) && !IsRootId(pid)) {
					if (IsSampleId(pid)) {
						if (!sampleToo) break;
					} else if (IsFolderId(pid)) {
						if (!folderToo) break;
					}
					c.add(pid);
					pid=(String)getParentId(pid);
				}
			}else {
				while (!found) {
					node=(ElementImpl )node.getParentNode();
					if (node==null)
						break;
					found=addGateOrSampleId(node, c, sampleToo);
				}
				if (found && sampleToo && folderToo) {
					pid=(String)getParentId(c.get(c.size()-1));
					while (!Basics.isEmpty(pid) && !ROOT_ID.equals(pid)) {
						c.add(pid);
						pid=(String)getParentId(c.get(c.size()-1));
					}
					final int N=c.size();
					for (int i=1;i<N;i++) 
						props.setProperty(c.get(i-1)+".parent", c.get(i));
				}
			}
			Collections.reverse(c);
		}
		return c;
	}
	
	public static String []GATE= {"Gate", "gating:*"};
	
	public String getGateType(final ElementImpl population) {
		if (population==null)
			return null;
		final ElementImpl node=GetFirstWild(population, GATE, 1);
		if (node==null)
			return null;
		return node.getNodeName().substring(7);
	}
	
	public Object[]getGateIdAndNodeAndType(final Object population) {
		if (!(population instanceof ElementImpl))
			return new Object[]{"", null, "?"};
		final ElementImpl node=GetFirstWild((ElementImpl)population, GATE, 1);
		if (node==null)
			return new Object[]{"", null, "?"};
		String id_=node.getAttribute("gating:id");
		if (Basics.isEmpty(id_)) 
			id_=((ElementImpl) node.getParentNode()).getAttribute("gating:id");
		return new Object[] {TYPE_GATE+":"+id_, node, node.getNodeName().substring(7)};
	}
	
	public Object getGateId(final ElementImpl population) {
		final ElementImpl node=GetFirstNode(population, "Gate");
		if (node==null)
			return null;
		return TYPE_GATE+":"+node.getAttribute("gating:id");
	}
	
	public Object getName(final Object nodeOrId) {
		final ElementImpl node;
		if (nodeOrId instanceof String) {
			node=getNodeById((String)nodeOrId);
		} else {
			node=(ElementImpl)nodeOrId;
		}
		if (node==null) return null;
		return node.getAttribute("name");
	}
	
	public boolean hasAnyChildren(final String pid) {
		final String prop=pid+".children";
		if (!props.containsKey(prop))
			getChildIdsString(pid); // cause XML items to be stored in props
		return !Basics.isEmpty(props.getProperty(prop));
	}

	public boolean rootIdsInSync=false;
	
	
	public Object getChildIdsString(final String pid){
		final String prop=pid+".children";
		final boolean justStarting=IsRootId(pid)&&!rootIdsInSync;
		if (!props.containsKey(prop) || justStarting) {
			// initialize props with children for pid
			props.setProperty(prop, "");
			if (IsSampleId(pid) || IsGateId(pid)) {//have gates
				final ElementImpl node=getNodeById(pid);
				if (node==null) return "";
				final Collection<ElementImpl> children=getSubPopulationList(node);
				for (final ElementImpl child:children) 
					setParent(getGateId(child), pid);
			} else if (justStarting)
				syncRootIds();
		}
		final String str=props.getProperty(prop);
		if (Basics.isEmpty(str)) return "";
		return str.trim();
	}
	
	public ArrayList<String>getSampleIdsInProperties(){
		return getSampleIdsInProperties(new ArrayList<String>(), ROOT_ID);
	}
	
	public ArrayList<String>getSampleIdsInProperties(ArrayList<String> c, final String pid){
		final String prop=pid+".children";
		if (props.containsKey(prop)) {
			final String str=this.props.getProperty(prop);
			if (!Basics.isEmpty(str)) {
				final String []ids=str.split(" ");
				for (final String id:ids) {
					if (IsSampleId(id)) {
						final int sampleNum=getSampleNumById(id);
						if (sampleNum<1)
							removeChild(pid, id);
						else
							c.add(id);
					}else 
						getSampleIdsInProperties(c,id);
				}
			}
		}
		return c;
	}

	public void syncRootIds() {
		final ArrayList<String>c=getSampleIdsInProperties();
		for (int i=1;i<=nSamples;i++) {
			final String sid=(String)getSampleIdByNum(i);
			if (!c.contains(sid))
				setParent(sid, ROOT_ID);
		}
		rootIdsInSync=true;
	}
	
	public Object getExpandableIds(final Collection<String>unexpandedIds) {
		final HashSet<String>pids=new HashSet<String>();
		final Collection<String>r=new ArrayList<String>();
		for (final String unexpandedId:unexpandedIds) {
			final String pid=(String)getParentId(unexpandedId);
			if (!pids.contains(pid)) {
				r.add(unexpandedId);
				pids.add(pid);
			}
		}
		final Iterator<String>it=r.iterator();
		while(it.hasNext()) {
			final String possibleAncestor=it.next();
			final Iterator<String> it2=r.iterator();
			while(it2.hasNext()) {
				final String possibleDescendent=it2.next();
				if (!possibleAncestor.equals(possibleDescendent)
						&& descendsFromAncestorOrItsSibling(
								possibleDescendent, possibleAncestor)) {
					it.remove();
					break;
				}
			}
		}
		String str=r.toString();
		str=str.substring(1, str.length()-1);
		return str.replace(", ", " ");
	}
	
	public boolean descendsFromAncestorOrItsSibling(
			final String suspectedDescendent, 
			final String suspectedAncestor) {
		String pid=(String)getParentId(suspectedDescendent);
		while(!Basics.isEmpty(pid)) {
			pid=(String)getParentId(pid);
			if (hasChild(pid, suspectedAncestor)) return true;
		}
		return false;
	}
	
	public boolean hasChild(final String pid, final String id) {
		if (Basics.isEmpty(pid))return false;
		String str=props.getProperty(pid+".children");
		if (str==null) {
			getChildIdsString(pid);
			str=props.getProperty(pid+".children");
		}
		return str.contains(id+" ");
	}
	
	
	public final static String [] Coords=new String [] {
			"gating:vertex", "gating:coordinate"};
	public static double [][]Reshape1to2(final double [] in){
		return Reshape1to2(in, null);
	}
	public static double [][]Reshape1to2(final double [] in, double [][]out){
		final int N=in.length;
		final int N2=N/2;
		if (out==null) {
			out=new double[N2][2];
		}
		int j=0;
		for(int i=0;i<N2;i++) {
			out[i][0]=in[j++];
			out[i][1]=in[j++];
		}
		return out;
	}
	
	public static double[][]GetPolygonCoords(final Object nodeObject){
		final double [] in=GetValueVector(nodeObject, Coords);
		if (in==null) return null;
		final int N=in.length;
		final int N2=N/2;
		final double [][]out;
		
		// full circle already?
		if (!(in[0]==in[N-2] && in[1]==in[N-1])) {
			// NO!!
			out=Reshape1to2(in, new double[N2+1][2]);
			Reshape1to2(in, out);
			out[N2][0]=out[0][0];
			out[N2][1]=out[0][1];
		}else {
			out=Reshape1to2(in);
		}
		return out;
	}

	public static final String [] Dims=new String[] {"gating:dimension"};
	public static final String []MinMax=new String [] {"gating:min", "gating:max"};
	public static double []GetRectCoords(final Object nodeObject){
		if (!(nodeObject instanceof ElementImpl)) return null;
		ElementImpl node=(ElementImpl)nodeObject;
		final List<ElementImpl>c=GetWildList(node, Dims);
		if (c.size()!=2)
			return null;
		final double []v=new double[4];
		int k=0;
		boolean isQuandrant=false;
		for (int j=0;j<2;j++) {
			for (int i=0;i<2;i++) {
				node=c.get(i);
				Double d;
				try {
					d=Double.valueOf(node.getAttribute(MinMax[j]));
				} catch (final Exception e) {
					d=null;
				}
				if (d==null) {
					d=Double.NaN;
					isQuandrant=true;
				}
				v[k]=d;
				k++;
			}
		}
		if (!isQuandrant) {
			v[2]=v[2]-v[0];
			v[3]=v[3]-v[1];
		}
		return v;
	}
	
	public final static String [] EllipseCoords=new String [] {
		"gating:edge", "gating:vertex", "gating:coordinate"};
	
	public static double[]GetEllipseCoords(final Object nodeObject){
		return GetValueVector(nodeObject, EllipseCoords);
	}
	
	public static double[]GetValueVector(final Object nodeObject, final String []names){
		return GetAttributeVector(nodeObject, names, "data-type:value");
	}
	
	public static double[]GetAttributeVector(
			final Object nodeObject, 
			final String []names, 
			final String attributeName){
		if (!(nodeObject instanceof ElementImpl)) return null;
		ElementImpl node=(ElementImpl)nodeObject;
		final Collection<ElementImpl>c=GetWildList(node, names);
		try {
			final int N=c.size();
			final double []v=new double[N];
			final Iterator<ElementImpl>it=c.iterator();
			for (int i=0;i<N;i++) {
				node=it.next();
				v[i]=Double.valueOf(node.getAttribute(attributeName));
			}
			return v;
		} catch( final Exception e ) {
			e.printStackTrace();
		}
		return null;
	}

	public static final String [] FcsDims=new String[] 
			{"gating:dimension", "data-type:fcs-dimension"};

	public static ElementImpl []GetFcsDims(final Object nodeObject){
		if (!(nodeObject instanceof ElementImpl)) return null;
		ElementImpl node=(ElementImpl)nodeObject;
		return GetWildList(node, FcsDims).toArray(emptyEle);
	}
	
	final static String []derivied=new String [] {"DerivedParameters", "DerivedParameter"};
	
	public Object getSampleDerived(final Object nodeOrNum){
		return getSampleItems(nodeOrNum, derivied, true);
	}
	
	public Object getSampleItems(final Object nodeOrNum, final String []names, final boolean all) {
		final ElementImpl in;
		if (nodeOrNum instanceof Integer) 
			in=(ElementImpl) sampleNodes.item(((Integer)nodeOrNum).intValue()-1);
		else if (nodeOrNum instanceof Double) 
			in=(ElementImpl) sampleNodes.item(((Double)nodeOrNum).intValue()-1);
		else if (nodeOrNum instanceof ElementImpl)
			in=(ElementImpl)nodeOrNum;
		else return null;
		if (all)
			return GetList(in, names);
		else
			return GetFirst(in, names, 1);
	}

	public int getSampleNumByName(final String name) {
		return getSampleNumByName(name, 0);
	}
	
	public int getSampleNumByName(final String name, final int defaultNum) {
		if (name.endsWith("*")) {
			final String name2=name.substring(0, name.length()-1);
			for (int i=1;i<=nSamples;i++) 
				if (getSampleNode(i).getAttribute("name").startsWith(name2))
					return i;
		} else {
			for (int i=1;i<=nSamples;i++) 
				if (getSampleNode(i).getAttribute("name").equals(name))
					return i;
		}
		return defaultNum;
	}

}
