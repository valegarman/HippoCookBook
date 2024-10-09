/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */

package edu.stanford.facs.swing;

/*
Java Swing, 2nd Edition
By Marc Loy, Robert Eckstein, Dave Wood, James Elliott, Brian Cole
ISBN: 0-596-00408-7
Publisher: O'Reilly 
 */
// TreeDragTest.java
//A simple starting point for testing the DnD JTree code.
//

import java.awt.Image;
import java.awt.Toolkit;
import java.awt.Point;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetContext;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

import javax.swing.AbstractButton;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import java.awt.datatransfer.Clipboard;


public class DndTrees extends AbstractButton {
	final List<File>missing=new ArrayList<File>();
	public static void ToClipBoard(final String file){
		DndTrees.ToClipBoard(file, false);
	}
	public static void ToClipBoard(final String file, final boolean asUrlFlavor){
		final Toolkit toolkit = Toolkit.getDefaultToolkit();
		final Clipboard clipboard = toolkit.getSystemClipboard();
		if (CpuInfo.isMac()){
			/* 
			 * there is a bug with java 7 and mac os x 10.7 and higher where
			 * only String data flavors can be copied to the clipboard not URL or File
			 * 
			 * The drag and drop of File or URL still works with java 7 and mac os.x.7.
			 * 
			 * To see proof that string copy to clipboard works you can set a breakpoint
			 * below.
			 */
			 StringSelection selection = new StringSelection("Java 7 clipboard copy works only for string");
			 clipboard.setContents(selection, selection);
		}
		if (asUrlFlavor){
			try{
				final String s="file://"+file;
				final URL url=new URL(s);
				final TransferableURL tf=new TransferableURL(url);
				clipboard.setContents(tf, tf);				
			} catch (final Exception e){
				e.printStackTrace();
			}
		}else{
			final ArrayList<File>f=new ArrayList<File>();
			final File fl=new File(file);
			final boolean exists=fl.exists();
			f.add(fl);
			final TransferableFile tf=new TransferableFile(f);
			clipboard.setContents(tf, tf);
		}
	}
	
	public String treeName="";
	public String noDropZone="<html><font color='red'>Drag over top of something valid...</font></html>";
	
	public DndTrees() {		
	}

	public static void reloadAll(final JTree tree){
		DefaultTreeModel model=(DefaultTreeModel)tree.getModel();
		final Enumeration e = tree.getExpandedDescendants(new TreePath (((DefaultMutableTreeNode)model.getRoot()).getPath()));
		model.reload();
		while (e != null && e.hasMoreElements()) {
			final TreePath tp=(TreePath)e.nextElement();
			tree.expandPath(tp); 
		}
	}
	
	public void manage(final JTree tree){
		manage(tree,null,5,5,10,10);
	}
	public void manage(final JTree tree, final Image dragImage, 
			final int dragOffSetX, final int dragOffSetY,
			final int dropOffSetX, final int dropOffSetY){
		new Drag(tree, DnDConstants.ACTION_COPY_OR_MOVE, dragImage, dragOffSetX, dragOffSetY);
		new Drop(tree, dropOffSetX, dropOffSetY);
	}
	
	public void rejectDragDrop(final String drag, final String drop){
		rejectDragDrop(drag, drop, false);
	}
	
	public void rejectDragDrop(final String drag, final String drop, final boolean ignoreCase){
		final DropRule dr=new DropRule(drag, drop, ignoreCase);
		dropRules.add(dr);
	}
	
	public void prohibitExportImport(final String prohibition){
		prohibitedExportImports.add(prohibition);
		
	}
	
	private List<String>prohibitedExportImports=new ArrayList<String>();
	private List<DropRule>dropRules=new ArrayList<DropRule>();
	
	static DefaultMutableTreeNode toTreeNode(final TreePath p){
		final Object c=p.getLastPathComponent();
		if (c instanceof DefaultMutableTreeNode){
			return (DefaultMutableTreeNode)c;
		}
		return null;
	}
	
	static String toTreeNodeUserObject(final TreeNode node){
		if (node instanceof DefaultMutableTreeNode){
			final Object c=((DefaultMutableTreeNode)node).getUserObject();
			if (c != null){
				return c.toString();
			}
		}
		return null;
	}

	
	class Drop implements DropTargetListener{
		private int dropOffSetX, dropOffSetY;
		private final JTree targetTree;

		Drop(final JTree tree, final int dropOffSetX, final int dropOffSetY) {
			targetTree = tree;
			this.dropOffSetX=dropOffSetX;
			this.dropOffSetY=dropOffSetY; 
			new DropTarget(targetTree, this);
		}

		private String stripAndFormat(final TreeNode node){
			if (node==null){
				return "";
			}
			String s=node.toString();
			final int htmlIdx=s.indexOf("<html>");
			if (htmlIdx>=0){
				s=s.substring(htmlIdx).replaceAll("<html>", "").replaceAll("</html>", "");
				//System.out.println("value="+s);
				return s;
			} else {
				final String cl=node.getClass().getName();
				final int clIdx=s.indexOf(cl);
				System.out.println("Class="+cl);
				System.out.println("clIdx="+clIdx+", toString()="+s);
				if (clIdx>=0 && clIdx<=s.length()-2){
					s=s.substring(clIdx+cl.length()+1);
				}
			}
			final int idx=s==null?-1:s.lastIndexOf("file://");
			if (idx>1 && s.charAt(idx-1)!='\''){
					s=s.substring(0, idx);
			}
			return s;
		}
		/*
		 * Drop Event Handlers
		 */
		private DefaultMutableTreeNode getNode(
				final DropTargetDragEvent dtde, final boolean showOverMsg) {
			final Point p = dtde.getLocation();
			final DropTargetContext dtc = dtde.getDropTargetContext();
			final JTree tree = (JTree) dtc.getComponent();
			final TreePath path = tree.getPathForLocation(p.x, p.y);
			final DefaultMutableTreeNode node;
			if (path==null){
				node=null;
			}else{
				node=(DefaultMutableTreeNode) path.getLastPathComponent();
			}
			if (showOverMsg && node != over){
				final String ss;
				final ToolTipOnDemand ttod=ToolTipOnDemand.getSingleton();
				if (node==null){
					ss=noDropZone;
				}else if ( node==lastDraggedNode){
					ss="<html><font color='red'>Can't drag & drop on self.</font></html>";
				}else{
					final String title;
					if (reject(node)){
						title="<font color='red'>You can't drop</font>";
					}else{
						title="You can drop";
					}
					final String different;
					if (tree != lastDraggedTree){
						different=" <i>(from different "+treeName+")</i> ";
					}else{
						different="";
					}
					lastDroppedText=stripAndFormat(node);
					lastDraggedText=stripAndFormat(lastDraggedNode);
					String s="<br>&nbsp;&nbsp;on to  \"<b>"+lastDroppedText+"</b>\"";
					ss="<html>" + title+
							"  \"<b>"+lastDraggedText+
							"</b>\""+different +s+"</html>";
					//System.out.println("Drop->"+lastDroppedText + "("+lastDroppedUserObjectString +")");
					//System.out.println("Drag->"+lastDraggedText + "("+lastDraggedUserObjectString+")");
					//System.out.println(ss);
				}
				ttod.close();
				ttod.show(tree, false,  p.x+dropOffSetX, p.y+dropOffSetY, null, ss);
				over=node;
			}
			return node;
		}
		
		private TreeNode over;
		
		boolean reject(final DefaultMutableTreeNode dropPath){
			if (dropRules.size()>0 && lastDraggedUserObjectString != null){
				final String dropUserObjectStr=toTreeNodeUserObject(dropPath);
				if (dropUserObjectStr != null){
					for (final DropRule dr:dropRules){
						if (dr.matches(lastDraggedUserObjectString, dropUserObjectStr)){
							return true;
						}
					}
				}
			}
			if (lastDraggedTree != this.targetTree){
				if (prohibitedExportImports.size()>0){
					final String dropUserObjectStr=toTreeNodeUserObject(dropPath);
					for (final String prohibition:prohibitedExportImports){
						if (lastDraggedNode==null){
							return true;
						}							
						if (lastDraggedUserObjectString.startsWith(prohibition)){
							return true;
						}
						if (dropUserObjectStr.startsWith(prohibition)){
							return true;
						}
					}
				}
			}
			return false;
		}


		public void dragEnter(final DropTargetDragEvent dtde) {
			final DefaultMutableTreeNode node = getNode(dtde,false);
			if (node==null || node==lastDraggedNode || reject(node)) {
				dtde.rejectDrag();
			} else {
				dtde.acceptDrag(dtde.getDropAction());
			}
		}

		public void dragOver(final DropTargetDragEvent dtde) {
			final DefaultMutableTreeNode node = getNode(dtde, true);
			if (node==null || node==lastDraggedNode || reject(node)) {
				dtde.rejectDrag();
			} else {
				dtde.acceptDrag(dtde.getDropAction());
			}
		}

		public void dragExit(DropTargetEvent dte) {
		}

		public void dropActionChanged(DropTargetDragEvent dtde) {
		}

		public void drop(final DropTargetDropEvent dtde) {
			final Point pt = dtde.getLocation();
			final DropTargetContext dtc = dtde.getDropTargetContext();
			final JTree tree = (JTree) dtc.getComponent();
			final TreePath parentPath = tree.getPathForLocation(pt.x, pt.y);
			if (parentPath==null){
				dtde.rejectDrop();
				return;
			}
			//System.out.println("Dropping on treeId="+tree.hashCode());
			final DefaultMutableTreeNode node = (DefaultMutableTreeNode) parentPath
					.getLastPathComponent();
			if (node==lastDraggedNode || reject(node)) {
				dtde.rejectDrop();
				return;
			}
			if (node==null ){
				System.out.println("ParentPath is NULL >> but we were over ... "+over);
				if (over ==null){
					return;
				}
				
			}
			try {
				final Transferable tr = dtde.getTransferable();
				final DataFlavor[] flavors = tr.getTransferDataFlavors();
				for (int i = 0; i < flavors.length; i++) {
					if (tr.isDataFlavorSupported(flavors[i])) {
						lastTransferable=tr;
						dtde.acceptDrop(dtde.getDropAction());
						lastDroppedPath=parentPath;
						lastDroppedTree=targetTree;
						lastDroppedNode=node;
						lastDroppedUserObjectString=toTreeNodeUserObject(lastDroppedNode);
						notifyDrop();
						dtde.dropComplete(true);
						//dbg.out("Drop occurred", dtde);
						return;
					}
				}
				dtde.rejectDrop();
			} catch (final Exception e) {
				e.printStackTrace();
				dtde.rejectDrop();
			}
		}
	}

	public static class DropData{
		public String key;
		public Transferable lastTransferable;
		public Object lastTransferableData;
		public JTree lastDraggedTree, lastDroppedTree;
		public TreePath lastDraggedPath, lastDroppedPath;
		public DefaultMutableTreeNode lastDraggedNode, lastDroppedNode;
		public String lastDraggedUserObjectString, lastDroppedUserObjectString;
	}
	public DropData lastDrop;
	public void notifyDrop() {
		lastDrop=new DropData();
		lastDrop.lastTransferable=lastTransferable;
		try {
			lastDrop.lastTransferableData=lastTransferable.getTransferData(DataFlavor.javaFileListFlavor);
		} catch (UnsupportedFlavorException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		lastDrop.lastDraggedNode=lastDraggedNode;
		lastDrop.lastDroppedNode=lastDroppedNode;
		lastDrop.lastDraggedTree=this.lastDraggedTree;
		lastDrop.lastDroppedTree=this.lastDroppedTree;
		lastDrop.lastDraggedPath=this.lastDraggedPath;
		lastDrop.lastDroppedPath=this.lastDroppedPath;
		lastDrop.lastDraggedUserObjectString=this.lastDraggedUserObjectString;
		lastDrop.lastDroppedUserObjectString=this.lastDroppedUserObjectString;
		System.out.println("lastDraggedNode="+lastDraggedNode);
		lastDrop.key="normal";
		final ActionEvent event=new ActionEvent(this, 22, "normal");
		fireActionPerformed(event);
		missing.clear();
	}

	public Transferable lastTransferable=null;
	public JTree lastDraggedTree, lastDroppedTree;
	public TreePath lastDraggedPath, lastDroppedPath;
	public DefaultMutableTreeNode lastDraggedNode, lastDroppedNode;
	public String lastDraggedUserObjectString, lastDroppedUserObjectString, lastDraggedText, lastDroppedText;
	
	class Drag implements DragSourceListener, DragGestureListener {
		private DragSource source;
		//private TransferableTreeNode transferable;
		private Transferable transferable;
		private JTree sourceTree;

		Drag(final JTree tree, final int actions, final Image dragImage,
				final int dragOffSetX, final int dragOffSetY) {
			sourceTree = tree;
			this.dragImage=dragImage;
			this.dragOffSetX=dragOffSetX;
			this.dragOffSetY=dragOffSetY;
			source = new DragSource();
			source.createDefaultDragGestureRecognizer(sourceTree,
					actions, this);
		}
		
		

		/*
		 * Drag Gesture Handler
		 */
		public void dragGestureRecognized(final DragGestureEvent dge) {
			final TreePath path = sourceTree.getSelectionPath();
			if ((path == null) || (path.getPathCount() <= 1)) {
				// We can't move the root node or an empty selection
				return;
			}
			lastDraggedTree=sourceTree;
			lastDraggedPath=path;
			lastDraggedNode=toTreeNode(path);
			lastDraggedUserObjectString=toTreeNodeUserObject(lastDraggedNode);
			lastDroppedPath=null;
			lastDroppedTree=null;
			lastDroppedNode=null;
			lastDroppedUserObjectString=null;
			//transferable = new TransferableTreeNode(lastDraggedUserObjectString);
			final List<File>files=new ArrayList<File>();
			final TreePath []tps=sourceTree.getSelectionPaths();
			for (final TreePath tp:tps){
				addFile(tp, files);
			}
			if (files.size()>0){
				System.out.println("Dragging "+files.size()+" files");
				transferable=new TransferableFile(files);
			}else{
				transferable=new TransferableTreeNode(lastDraggedUserObjectString);
			}
			if (dragImage==null){
				source.startDrag(dge, null, transferable, this);
			}else{
				source.startDrag(dge, null, dragImage, new Point(dragOffSetX, dragOffSetY), transferable, this);
			}
			missing.clear();
			for (final File f:files){
				if (!f.exists()){
					missing.add(f);
				}
			}
		}

		private Image dragImage=null;
		private int dragOffSetX, dragOffSetY;

		/*
		 * Drag Event Handlers
		 */
		public void dragEnter(DragSourceDragEvent dsde) {
		}

		public void dragExit(DragSourceEvent dse) {
		}

		public void dragOver(DragSourceDragEvent dsde) {
		}

		public void dropActionChanged(DragSourceDragEvent dsde) {
			System.out.println("Action: " + dsde.getDropAction());
			System.out.println("Target Action: " + dsde.getTargetActions());
			System.out.println("User Action: " + dsde.getUserAction());
		}

		public void dragDropEnd(final DragSourceDropEvent dsde) {
			System.out.println("Drop Action: " + dsde.getDropAction());
			ToolTipOnDemand.getSingleton().close();
			if (missing.size()>0){
				lastDrop=new DropData();
				lastDrop.lastDraggedNode=lastDraggedNode;
				lastDrop.lastDroppedNode=lastDroppedNode;
				lastDrop.lastDraggedTree=lastDraggedTree;
				lastDrop.lastDroppedTree=lastDroppedTree;
				lastDrop.lastDraggedPath=lastDraggedPath;
				lastDrop.lastDroppedPath=lastDroppedPath;
				lastDrop.lastDraggedUserObjectString=lastDraggedUserObjectString;
				lastDrop.lastDroppedUserObjectString=lastDroppedUserObjectString;
				lastDrop.key="Missing";
				final ActionEvent event=new ActionEvent(this, 22, "missing "+missing.iterator().next().getAbsolutePath());
				fireActionPerformed(event);
			}
			lastDraggedNode=null;
			
		}
		
		boolean addFile(final Object input, final List<File>files){
			boolean added=false;
			final String str;
			if (input instanceof String){
				str=(String)input;
			}else{
				str=toTreeNodeUserObject(toTreeNode((TreePath)input));
			}
			if (str!=null){
				final int idx=str==null?-1:str.lastIndexOf("file://");
				if (idx>=0){
					final String fp=str.substring(idx+7);
					System.out.println(fp);
					files.add(new File(fp));
					added=true;
				} 
			}
			return added;
		}
	}

	static class TransferableTreeNode implements Transferable {
		private static DataFlavor TREE_PATH_FLAVOR = new DataFlavor(String.class,
			      "Tree Path");
		private DataFlavor flavors[] = { TREE_PATH_FLAVOR };
		private String path;

		public TransferableTreeNode(String tp) {
			path = tp;
		}

		public synchronized DataFlavor[] getTransferDataFlavors() {
			return flavors;
		}

		public boolean isDataFlavorSupported(final DataFlavor flavor) {
			return (flavor.getRepresentationClass() == String.class);
		}

		public synchronized Object getTransferData(DataFlavor flavor)
				throws UnsupportedFlavorException, IOException {
			if (isDataFlavorSupported(flavor)) {
				return (Object) path;
			} else {
				throw new UnsupportedFlavorException(flavor);
			}
		}
	}
	
	static class TransferableFile implements Transferable, ClipboardOwner {
		private final DataFlavor flavors[] = { DataFlavor.javaFileListFlavor };
		private final List<File> files;

	    public void lostOwnership(Clipboard clipboard, Transferable contents) {
	    }

		public TransferableFile(final List<File>files) {
			this.files=files;
		}

		public synchronized DataFlavor[] getTransferDataFlavors() {
			return flavors;
		}

		public boolean isDataFlavorSupported(final DataFlavor flavor) {
			return (flavor.getRepresentationClass() == java.util.List.class);
		}

		public synchronized Object getTransferData(final DataFlavor flavor)
				throws UnsupportedFlavorException, IOException {
			if (isDataFlavorSupported(flavor)) {
				return (Object) files;
			} else {
				throw new UnsupportedFlavorException(flavor);
			}
		}
	}

	static class TransferableURL implements Transferable, ClipboardOwner {

		private final DataFlavor urlFlavor =
				new DataFlavor("application/x-java-url; class=java.net.URL",
				"Image URL");
		private final DataFlavor flavors[] = { urlFlavor };
		private final java.net.URL url;

	    public void lostOwnership(Clipboard clipboard, Transferable contents) {
	    }

		public TransferableURL(final java.net.URL url) {
			this.url=url;
		}

		public synchronized DataFlavor[] getTransferDataFlavors() {
			return flavors;
		}

		public boolean isDataFlavorSupported(final DataFlavor flavor) {
			return (flavor.getRepresentationClass() == java.net.URL.class);
		}
		public synchronized Object getTransferData(final DataFlavor flavor)
				throws UnsupportedFlavorException, IOException {
			if (isDataFlavorSupported(flavor)) {
				return (Object) url;
			} else {
				throw new UnsupportedFlavorException(flavor);
			}
		}
	}

}

