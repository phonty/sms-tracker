/*
 * Copyright (C) 2018 David Barry <david.barry at crick dot ac dot uk>
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

import java.awt.Component;
import java.awt.Container;
import java.util.Properties;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.JToggleButton;

/**
 *
 * @author David Barry <david.barry at crick dot ac dot uk>
 */
public class PropertyExtractor {

    public static void setProperties(Properties props, Container container) {
        Component[] comps = container.getComponents();
        if (props == null) {
            props = new Properties();
        }
        for (Component c : comps) {
            if (c instanceof Container) {
                setProperties(props, (Container) c);
            }
            if (c instanceof JLabel) {
                if (((JLabel) c).getLabelFor() instanceof JTextField) {
                    props.setProperty(((JLabel) c).getText(), ((JTextField) ((JLabel) c).getLabelFor()).getText());
                } else if (((JLabel) c).getLabelFor() instanceof JComboBox) {
                    props.setProperty(((JLabel) c).getText(), ((JComboBox) ((JLabel) c).getLabelFor()).getSelectedItem().toString());
                }
            } else if (c instanceof JToggleButton) {
                props.setProperty(((JToggleButton) c).getText(), String.format("%b", ((JToggleButton) c).isSelected()));
            }
        }
    }
}
