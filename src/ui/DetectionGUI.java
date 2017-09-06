/* 
 * Copyright (C) 2014 David Barry <david.barry at cancer.org.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package ui;

import IAClasses.Utils;
import Particle_Analysis.Particle_Tracker;
import ParticleTracking.GPUAnalyse;
import Particle.Particle;
import Particle.ParticleArray;
import ParticleTracking.UserVariables;
import UIClasses.UIMethods;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.process.ImageProcessor;
import java.awt.Canvas;
import java.awt.Color;
import java.awt.Component;
import java.util.ArrayList;
import java.util.Properties;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JToggleButton;

public class DetectionGUI extends javax.swing.JDialog {

    protected final Particle_Tracker analyser;
    protected final ImagePlus imp;
    protected final String title;
    protected boolean wasOKed = false, monoChrome;
    protected static final String spatResLabelText = "Spatial resolution (" + IJ.micronSymbol + "m/pixel):";
    protected static final String chan1MaxThreshLabelText = "C1 minimum peak size:";
    protected static final String chan2MaxThreshLabelText = "C2 minimum peak size:";
    protected static final String c1CurveFitTolLabelText = "Curve fit tolerance:";
    protected static final String preprocessToggleText = "Pre-Process Images";
    protected static final String gpuToggleText = "Use GPU";
    protected static final String redSigEstText = "PSF radius (" + IJ.micronSymbol + "m):";
    protected static final String greenSigEstText = "C2 PSF Width (" + IJ.micronSymbol + "m):";
    protected static final String DETECT_MODE = "Detection Mode:";
    protected static final String blobSizeText = "Blob size (" + IJ.micronSymbol + "m):";
    protected static final String filterRadiusText = "Gaussian Filter Radius (" + IJ.micronSymbol + "m):";
    protected static final DefaultComboBoxModel<String> DETECT_MODE_OPTIONS = new DefaultComboBoxModel(new String[]{"Points", "Blobs", "PSFs"});
    private Properties props;

    /**
     * Creates new form UserInterface
     */
    public DetectionGUI(java.awt.Frame parent, boolean modal, String title, Particle_Tracker analyser, boolean monoChrome) {
        super(parent, modal);
        this.title = title;
        this.analyser = analyser;
        ImageStack[] stacks = analyser.getStacks();
        this.monoChrome = monoChrome;
        if (this.monoChrome) {
            stacks[1] = null;
        }
        this.imp = new ImagePlus("", Utils.updateImage(stacks[0], stacks[1], 1));
        if (monoChrome) {
            UserVariables.setColocal(!monoChrome);
        }
        initComponents();
        UIMethods.centreDialog(this);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {
        java.awt.GridBagConstraints gridBagConstraints;

        jPanel3 = new javax.swing.JPanel();
        okButton = new javax.swing.JButton();
        cancelButton = new javax.swing.JButton();
        detectionPanel = new javax.swing.JPanel();
        spatResLabel = new javax.swing.JLabel();
        chan1MaxThreshLabel = new javax.swing.JLabel();
        spatResTextField = new javax.swing.JTextField();
        chan1MaxThreshTextField = new javax.swing.JTextField();
        preprocessToggleButton = new javax.swing.JToggleButton();
        curveFitTolLabel = new javax.swing.JLabel();
        curveFitTolTextField = new javax.swing.JTextField();
        gpuToggleButton = new javax.swing.JToggleButton();
        sigmaLabel = new javax.swing.JLabel();
        sigmaTextField = new javax.swing.JTextField();
        chan2MaxThreshLabel = new javax.swing.JLabel();
        chan2MaxThreshTextField = new javax.swing.JTextField();
        detectionModeComboBox = new javax.swing.JComboBox<>();
        detectionModeLabel = new javax.swing.JLabel();
        blobSizeLabel = new javax.swing.JLabel();
        blobSizeTextField = new javax.swing.JTextField();
        filterRadiusLabel = new javax.swing.JLabel();
        filterRadiusTextField = new javax.swing.JTextField();
        jPanel2 = new javax.swing.JPanel();
        canvas1 = new ImageCanvas(imp);
        previewTextField = new javax.swing.JTextField();
        previewToggleButton = new javax.swing.JToggleButton();
        previewScrollBar = new java.awt.Scrollbar();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
        setTitle(title);
        getContentPane().setLayout(new java.awt.GridBagLayout());

        jPanel3.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jPanel3.setLayout(new java.awt.GridBagLayout());

        okButton.setText("Run");
        okButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                okButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
        jPanel3.add(okButton, gridBagConstraints);

        cancelButton.setText("Cancel");
        cancelButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                cancelButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
        jPanel3.add(cancelButton, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.2;
        getContentPane().add(jPanel3, gridBagConstraints);

        detectionPanel.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        detectionPanel.setLayout(new java.awt.GridBagLayout());

        spatResLabel.setText(spatResLabelText);
        spatResLabel.setLabelFor(spatResTextField);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(spatResLabel, gridBagConstraints);

        chan1MaxThreshLabel.setText(chan1MaxThreshLabelText);
        chan1MaxThreshLabel.setLabelFor(chan1MaxThreshTextField);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(chan1MaxThreshLabel, gridBagConstraints);

        spatResTextField.setText(String.valueOf(UserVariables.getSpatialRes()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(spatResTextField, gridBagConstraints);

        chan1MaxThreshTextField.setText(String.valueOf(UserVariables.getChan1MaxThresh()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(chan1MaxThreshTextField, gridBagConstraints);

        preprocessToggleButton.setText(preprocessToggleText);
        preprocessToggleButton.setSelected(UserVariables.isPreProcess());
        preprocessToggleButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                preprocessToggleButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 10;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(preprocessToggleButton, gridBagConstraints);

        curveFitTolLabel.setText(c1CurveFitTolLabelText);
        curveFitTolLabel.setLabelFor(curveFitTolTextField);
        curveFitTolLabel.setEnabled(UserVariables.getDetectionMode()==UserVariables.GAUSS);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 9;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(curveFitTolLabel, gridBagConstraints);

        curveFitTolTextField.setText(String.valueOf(UserVariables.getCurveFitTol()));
        curveFitTolTextField.setEnabled(UserVariables.getDetectionMode()==UserVariables.GAUSS);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 9;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(curveFitTolTextField, gridBagConstraints);

        gpuToggleButton.setText(gpuToggleText);
        gpuToggleButton.setEnabled(analyser.isGpuEnabled());
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 12;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(gpuToggleButton, gridBagConstraints);

        sigmaLabel.setText(redSigEstText);
        sigmaLabel.setLabelFor(sigmaTextField);
        sigmaLabel.setEnabled(UserVariables.getDetectionMode()==UserVariables.GAUSS);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 8;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(sigmaLabel, gridBagConstraints);

        sigmaTextField.setText(String.valueOf(UserVariables.getSigEstRed()));
        sigmaTextField.setEnabled(UserVariables.getDetectionMode()==UserVariables.GAUSS);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 8;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(sigmaTextField, gridBagConstraints);

        chan2MaxThreshLabel.setText(chan2MaxThreshLabelText);
        chan2MaxThreshLabel.setLabelFor(chan2MaxThreshTextField);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(chan2MaxThreshLabel, gridBagConstraints);

        chan2MaxThreshTextField.setText(String.valueOf(UserVariables.getChan2MaxThresh()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 3;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(chan2MaxThreshTextField, gridBagConstraints);

        detectionModeComboBox.setModel(DETECT_MODE_OPTIONS);
        detectionModeComboBox.setSelectedIndex(UserVariables.getDetectionMode()-UserVariables.MAXIMA);
        detectionModeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                detectionModeComboBoxActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(detectionModeComboBox, gridBagConstraints);

        detectionModeLabel.setText(DETECT_MODE);
        detectionModeLabel.setLabelFor(detectionModeComboBox);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(detectionModeLabel, gridBagConstraints);

        blobSizeLabel.setText(blobSizeText);
        blobSizeLabel.setLabelFor(blobSizeTextField);
        blobSizeLabel.setEnabled(!(UserVariables.getDetectionMode()==UserVariables.GAUSS));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(blobSizeLabel, gridBagConstraints);

        blobSizeTextField.setText(String.format("%1.3f", UserVariables.getBlobSize()));
        blobSizeTextField.setEnabled(!(UserVariables.getDetectionMode()==UserVariables.GAUSS));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 6;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 1.0;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(blobSizeTextField, gridBagConstraints);

        filterRadiusLabel.setText(filterRadiusText);
        filterRadiusLabel.setLabelFor(filterRadiusTextField);
        filterRadiusLabel.setEnabled(UserVariables.isPreProcess() && !(UserVariables.getDetectionMode()==UserVariables.BLOBS));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 0;
        gridBagConstraints.gridy = 11;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 0, 10);
        detectionPanel.add(filterRadiusLabel, gridBagConstraints);

        filterRadiusTextField.setText(String.format("%1.3f", UserVariables.getFilterRadius()));
        filterRadiusTextField.setEnabled(UserVariables.isPreProcess()&&!(UserVariables.getDetectionMode()==UserVariables.BLOBS));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 11;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_END;
        gridBagConstraints.insets = new java.awt.Insets(0, 0, 0, 10);
        detectionPanel.add(filterRadiusTextField, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 1.0;
        getContentPane().add(detectionPanel, gridBagConstraints);

        jPanel2.setBorder(javax.swing.BorderFactory.createEtchedBorder());
        jPanel2.setLayout(new java.awt.GridBagLayout());

        canvas1.setMinimumSize(new java.awt.Dimension(analyser.getStacks()[0].getWidth()/4,analyser.getStacks()[0].getHeight()/4));
        canvas1.setPreferredSize(new java.awt.Dimension(analyser.getStacks()[0].getWidth(),analyser.getStacks()[0].getHeight()));
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 0;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.8;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 0, 10);
        jPanel2.add(canvas1, gridBagConstraints);

        previewTextField.setText(String.valueOf(previewScrollBar.getValue()));
        previewTextField.setEditable(false);
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 2;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
        gridBagConstraints.weightx = 0.2;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 10, 10);
        jPanel2.add(previewTextField, gridBagConstraints);

        previewToggleButton.setText("Preview");
        previewToggleButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                previewToggleButtonActionPerformed(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 1;
        gridBagConstraints.gridwidth = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.weightx = 1.0;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(10, 10, 10, 10);
        jPanel2.add(previewToggleButton, gridBagConstraints);

        previewScrollBar.setOrientation(java.awt.Scrollbar.HORIZONTAL);
        previewScrollBar.setValues(1, 1, 1, analyser.getStacks()[0].getSize());
        previewScrollBar.addAdjustmentListener(new java.awt.event.AdjustmentListener() {
            public void adjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {
                previewScrollBarAdjustmentValueChanged(evt);
            }
        });
        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.gridx = 1;
        gridBagConstraints.gridy = 2;
        gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
        gridBagConstraints.anchor = java.awt.GridBagConstraints.WEST;
        gridBagConstraints.weightx = 0.8;
        gridBagConstraints.weighty = 0.1;
        gridBagConstraints.insets = new java.awt.Insets(0, 10, 10, 10);
        jPanel2.add(previewScrollBar, gridBagConstraints);

        gridBagConstraints = new java.awt.GridBagConstraints();
        gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
        gridBagConstraints.weightx = 0.5;
        gridBagConstraints.weighty = 1.0;
        getContentPane().add(jPanel2, gridBagConstraints);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void previewToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_previewToggleButtonActionPerformed
        previewScrollBarAdjustmentValueChanged(null);
    }//GEN-LAST:event_previewToggleButtonActionPerformed

    private void cancelButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_cancelButtonActionPerformed
        this.dispose();
        wasOKed = false;
    }//GEN-LAST:event_cancelButtonActionPerformed

    private void okButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_okButtonActionPerformed
        if (!setVariables()) {
            return;
        }
        this.dispose();
        wasOKed = true;
    }//GEN-LAST:event_okButtonActionPerformed

    private void previewScrollBarAdjustmentValueChanged(java.awt.event.AdjustmentEvent evt) {//GEN-FIRST:event_previewScrollBarAdjustmentValueChanged
        previewTextField.setText(String.valueOf(previewScrollBar.getValue()));
        if (previewToggleButton.isSelected() && !previewScrollBar.getValueIsAdjusting() && setVariables()) {
            viewDetections(analyser, monoChrome, Double.parseDouble(spatResTextField.getText()), previewScrollBar.getValue(), canvas1, imp);
        }
    }//GEN-LAST:event_previewScrollBarAdjustmentValueChanged

    private void detectionModeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_detectionModeComboBoxActionPerformed
        setVariables();
        boolean psfs = UserVariables.getDetectionMode() == UserVariables.GAUSS;
        curveFitTolTextField.setEnabled(psfs);
        curveFitTolLabel.setEnabled(psfs);
        sigmaTextField.setEnabled(psfs);
        sigmaLabel.setEnabled(psfs);
        blobSizeLabel.setEnabled(!psfs);
        blobSizeTextField.setEnabled(!psfs);
        preprocessToggleButtonActionPerformed(evt);
    }//GEN-LAST:event_detectionModeComboBoxActionPerformed

    private void preprocessToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_preprocessToggleButtonActionPerformed
        preprocessToggleButton.setEnabled(!(UserVariables.getDetectionMode() == UserVariables.BLOBS));
        filterRadiusLabel.setEnabled(preprocessToggleButton.isSelected() && !(UserVariables.getDetectionMode() == UserVariables.BLOBS));
        filterRadiusTextField.setEnabled(preprocessToggleButton.isSelected() && !(UserVariables.getDetectionMode() == UserVariables.BLOBS));
    }//GEN-LAST:event_preprocessToggleButtonActionPerformed

    boolean setVariables() {
        try {
            UserVariables.setDetectionMode(detectionModeComboBox.getSelectedIndex() + UserVariables.MAXIMA);
            UserVariables.setChan1MaxThresh(Double.parseDouble(chan1MaxThreshTextField.getText()));
            UserVariables.setChan2MaxThresh(Double.parseDouble(chan2MaxThreshTextField.getText()));
            UserVariables.setSpatialRes(Double.parseDouble(spatResTextField.getText()));
            UserVariables.setCurveFitTol(Double.parseDouble(curveFitTolTextField.getText()));
            UserVariables.setPreProcess(preprocessToggleButton.isSelected());
            UserVariables.setGpu(gpuToggleButton.isSelected());
            UserVariables.setSigEstRed(Double.parseDouble(sigmaTextField.getText()));
//            UserVariables.setSigEstGreen(Double.parseDouble(greenSigmaTextField.getText()));
            UserVariables.setBlobSize(Double.parseDouble(blobSizeTextField.getText()));
            UserVariables.setFilterRadius(Double.parseDouble(filterRadiusTextField.getText()));
        } catch (NumberFormatException e) {
            IJ.error("Number formatting error " + e.toString());
            return false;
        }
        setProperties();
        return true;
    }

    private void setProperties() {
        Component[] comps = detectionPanel.getComponents();
        props = new Properties();
        for (Component c : comps) {
            if (c instanceof JLabel && ((JLabel) c).getLabelFor() instanceof JTextField) {
                props.setProperty(((JLabel) c).getText(), ((JTextField) ((JLabel) c).getLabelFor()).getText());
            } else if (c instanceof JToggleButton) {
                props.setProperty(((JToggleButton) c).getText(), String.format("%b", ((JToggleButton) c).isSelected()));
            }
        }
    }

    public static void viewDetections(Particle_Tracker analyser, boolean monoChrome, double spatRes, int psv, Canvas canvas1, ImagePlus imp) {
        analyser.calcParticleRadius(UserVariables.getSpatialRes());
        ImageStack stacks[] = analyser.getStacks();
        if (monoChrome) {
            stacks[1] = null;
        }
        ParticleArray detections;
        if (psv < 1) {
            psv = 1;
        }
        if (analyser instanceof GPUAnalyse && UserVariables.isGpu()) {
            detections = ((GPUAnalyse) analyser).cudaFindParticles(false, psv - 1, psv - 1, stacks[1]);
        } else {
            detections = analyser.findParticles(false, psv - 1, psv - 1, UserVariables.getCurveFitTol(), stacks[0], stacks[1]);
        }
        if (detections != null) {
            ImageProcessor output = Utils.updateImage(stacks[0], stacks[1], psv);
            double mag = 1.0 / UIMethods.getMagnification(output, canvas1);
//            double sr = 1.0 / spatRes;
//            int radius = analyser.calcParticleRadius(UserVariables.getSpatialRes());
            ArrayList<Particle> particles = detections.getLevel(0);
            Color c1Color = !monoChrome ? Color.red : Color.white;
            Color c2Color = !monoChrome ? Color.green : Color.white;
            output.setLineWidth(1);
            for (Particle p1 : particles) {
                Particle p2 = p1.getColocalisedParticle();
                drawParticle(output, true, c1Color, p1);
                if (p2 != null) {
                    drawParticle(output, true, c2Color, p2);
                }
            }
            imp.setProcessor("", output);
            ((ImageCanvas) canvas1).setMagnification(mag);
            canvas1.repaint();
        }
    }

    public static void drawParticle(ImageProcessor image,
            boolean drawOval, Color colour, Particle p) {
        image.setColor(colour);
        int radius = (int) Math.round(UserVariables.getBlobSize() / UserVariables.getSpatialRes());
        int x = (int) Math.round(p.getX() / UserVariables.getSpatialRes());
        int y = (int) Math.round(p.getY() / UserVariables.getSpatialRes());
        switch (UserVariables.getDetectionMode()) {
            case UserVariables.BLOBS:
                image.drawOval((x - radius), (y - radius), 2 * radius, 2 * radius);
                break;
            case UserVariables.GAUSS:
                radius = (int) Math.round(UserVariables.getSigEstRed() / UserVariables.getSpatialRes());
                image.drawOval((x - radius), (y - radius), 2 * radius, 2 * radius);
                break;
            default:
                image.drawLine(x, y - radius, x, y + radius);
                image.drawLine(x + radius, y, x - radius, y);
        }
    }

    public boolean isWasOKed() {
        return wasOKed;
    }

    public static String getSpatResLabelText() {
        return spatResLabelText;
    }

    public static String getChan1MaxThreshLabelText() {
        return chan1MaxThreshLabelText;
    }

    public static String getCurveFitTolLabelText() {
        return c1CurveFitTolLabelText;
    }

    public static String getPreprocessToggleText() {
        return preprocessToggleText;
    }

    public static String getGpuToggleText() {
        return gpuToggleText;
    }

    public static String getRedSigEstText() {
        return redSigEstText;
    }

    public static String getGreenSigEstText() {
        return greenSigEstText;
    }

    public static String getChan2MaxThreshLabelText() {
        return chan2MaxThreshLabelText;
    }

    public Properties getProperties() {
        return props;
    }


    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel blobSizeLabel;
    private javax.swing.JTextField blobSizeTextField;
    private javax.swing.JButton cancelButton;
    private java.awt.Canvas canvas1;
    private javax.swing.JLabel chan1MaxThreshLabel;
    private javax.swing.JTextField chan1MaxThreshTextField;
    private javax.swing.JLabel chan2MaxThreshLabel;
    private javax.swing.JTextField chan2MaxThreshTextField;
    private javax.swing.JLabel curveFitTolLabel;
    private javax.swing.JTextField curveFitTolTextField;
    private javax.swing.JComboBox<String> detectionModeComboBox;
    private javax.swing.JLabel detectionModeLabel;
    private javax.swing.JPanel detectionPanel;
    private javax.swing.JLabel filterRadiusLabel;
    private javax.swing.JTextField filterRadiusTextField;
    private javax.swing.JToggleButton gpuToggleButton;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JButton okButton;
    private javax.swing.JToggleButton preprocessToggleButton;
    private java.awt.Scrollbar previewScrollBar;
    private javax.swing.JTextField previewTextField;
    private javax.swing.JToggleButton previewToggleButton;
    private javax.swing.JLabel sigmaLabel;
    private javax.swing.JTextField sigmaTextField;
    private javax.swing.JLabel spatResLabel;
    private javax.swing.JTextField spatResTextField;
    // End of variables declaration//GEN-END:variables
}
