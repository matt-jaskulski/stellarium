/*
 * Stellarium Satellites Plug-in: satellites custom filter feature
 * Copyright (C) 2022 Alexander Wolf
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include "SatellitesFilterDialog.hpp"
#include "ui_satellitesFilterDialog.h"

#include "StelApp.hpp"
#include "StelTranslator.hpp"

SatellitesFilterDialog::SatellitesFilterDialog()
	: StelDialog("SatellitesFilter")
{
	ui = new Ui_satellitesFilterDialog;
}

SatellitesFilterDialog::~SatellitesFilterDialog()
{
	delete ui;
}

void SatellitesFilterDialog::retranslate()
{
	if (dialog)
	{
		ui->retranslateUi(dialog);
		populateTexts();
	}
}

void SatellitesFilterDialog::setVisible(bool visible)
{
	StelDialog::setVisible(visible);
}

void SatellitesFilterDialog::createDialogContent()
{
	ui->setupUi(dialog);

	connect(ui->closeStelWindow, SIGNAL(clicked()), this, SLOT(close()));
	connect(ui->TitleBar, SIGNAL(movedTo(QPoint)), this, SLOT(handleMovedTo(QPoint)));

	connectBoolProperty(ui->inclinationCheckBox,  "Satellites.flagCFInclination");
	connectDoubleProperty(ui->minInclination,     "Satellites.minCFInclination");
	connectDoubleProperty(ui->maxInclination,     "Satellites.maxCFInclination");
	connectBoolProperty(ui->periodCheckBox,       "Satellites.flagCFPeriod");
	connectDoubleProperty(ui->minPeriod,          "Satellites.minCFPeriod");
	connectDoubleProperty(ui->maxPeriod,          "Satellites.maxCFPeriod");
	connectBoolProperty(ui->eccentricityCheckBox, "Satellites.flagCFEccentricity");
	connectDoubleProperty(ui->minEccentricity,    "Satellites.minCFEccentricity");
	connectDoubleProperty(ui->maxEccentricity,    "Satellites.maxCFEccentricity");
	connectBoolProperty(ui->apogeeCheckBox,       "Satellites.flagCFApogee");
	connectDoubleProperty(ui->minApogee,          "Satellites.minCFApogee");
	connectDoubleProperty(ui->maxApogee,          "Satellites.maxCFApogee");
	connectBoolProperty(ui->perigeeCheckBox,      "Satellites.flagCFPerigee");
	connectDoubleProperty(ui->minPerigee,         "Satellites.minCFPerigee");
	connectDoubleProperty(ui->maxPerigee,         "Satellites.maxCFPerigee");
	connectBoolProperty(ui->altitudeCheckBox,     "Satellites.flagCFAltitude");
	connectDoubleProperty(ui->minAltitude,        "Satellites.minCFAltitude");
	connectDoubleProperty(ui->maxAltitude,        "Satellites.maxCFAltitude");
	connectBoolProperty(ui->rcsCheckBox,          "Satellites.flagCFRCS");
	connectDoubleProperty(ui->minRCS,             "Satellites.minCFRCS");
	connectDoubleProperty(ui->maxRCS,             "Satellites.maxCFRCS");
	connectBoolProperty(ui->stdMagnitudeCheckBox, "Satellites.flagCFKnownStdMagnitude");

	updateMinMaxInclination(ui->inclinationCheckBox->isChecked());
	connect(ui->inclinationCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxInclination(bool)));
	updateMinMaxApogee(ui->apogeeCheckBox->isChecked());
	connect(ui->apogeeCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxApogee(bool)));
	updateMinMaxPerigee(ui->perigeeCheckBox->isChecked());
	connect(ui->perigeeCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxPerigee(bool)));
	updateMinMaxAltitude(ui->altitudeCheckBox->isChecked());
	connect(ui->altitudeCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxAltitude(bool)));
	updateMinMaxPeriod(ui->periodCheckBox->isChecked());
	connect(ui->periodCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxPeriod(bool)));
	updateMinMaxEccentricity(ui->eccentricityCheckBox->isChecked());
	connect(ui->eccentricityCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxEccentricity(bool)));
	updateMinMaxRCS(ui->rcsCheckBox->isChecked());
	connect(ui->rcsCheckBox, SIGNAL(clicked(bool)), this, SLOT(updateMinMaxRCS(bool)));

	populateTexts();
}

void SatellitesFilterDialog::updateMinMaxInclination(bool state)
{
	ui->minInclination->setEnabled(state);
	ui->maxInclination->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxApogee(bool state)
{
	ui->minApogee->setEnabled(state);
	ui->maxApogee->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxPerigee(bool state)
{
	ui->minPerigee->setEnabled(state);
	ui->maxPerigee->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxAltitude(bool state)
{
	ui->minAltitude->setEnabled(state);
	ui->maxAltitude->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxPeriod(bool state)
{
	ui->minPeriod->setEnabled(state);
	ui->maxPeriod->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxEccentricity(bool state)
{
	ui->minEccentricity->setEnabled(state);
	ui->maxEccentricity->setEnabled(state);
}

void SatellitesFilterDialog::updateMinMaxRCS(bool state)
{
	ui->minRCS->setEnabled(state);
	ui->maxRCS->setEnabled(state);
}

void SatellitesFilterDialog::populateTexts()
{
	QString km = qc_("km", "distance");
	QString m = qc_("m", "distance");
	QString min = qc_("m", "time");
	ui->minApogee->setSuffix(QString(" %1").arg(km));
	ui->maxApogee->setSuffix(QString(" %1").arg(km));
	ui->minPerigee->setSuffix(QString(" %1").arg(km));
	ui->maxPerigee->setSuffix(QString(" %1").arg(km));
	ui->altitudeCheckBox->setToolTip(q_("The satellite can be located in selected range of altitudes"));
	ui->minAltitude->setSuffix(QString(" %1").arg(km));
	ui->maxAltitude->setSuffix(QString(" %1").arg(km));
	ui->minPeriod->setSuffix(QString(" %1").arg(min));
	ui->maxPeriod->setSuffix(QString(" %1").arg(min));
	ui->minInclination->setSuffix("°");
	ui->maxInclination->setSuffix("°");
	ui->rcsCheckBox->setToolTip(q_("Radar cross-section"));
	ui->minRCS->setSuffix(QString(" %1%2").arg(m, "²"));
	ui->maxRCS->setSuffix(QString(" %1%2").arg(m, "²"));
}
