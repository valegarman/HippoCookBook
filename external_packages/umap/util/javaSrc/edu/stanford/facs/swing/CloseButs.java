/***
 * Author: Stephen Meehan, swmeehan@stanford.edu
 * 
 * Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
 * 
 * License: BSD 3 clause
 * 
 */

package edu.stanford.facs.swing;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JPanel;

public class CloseButs {
	public boolean cancelled=true;
	public final JButton done=new JButton("Ok"),cancel=new JButton("Cancel");
    public CloseButs(final JDialog dlg, final JPanel buts){
		final JButton cancel=new JButton("Cancel");
		CpuInfo.registerEscape(dlg, cancel);
		buts.add(cancel);
		buts.add(done);
		done.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cancelled=false;
				dlg.dispose();
			}
		});
		cancel.addActionListener(CpuInfo.getCloseAction(dlg));
		
		dlg.getRootPane().setDefaultButton(done);
    }
}
