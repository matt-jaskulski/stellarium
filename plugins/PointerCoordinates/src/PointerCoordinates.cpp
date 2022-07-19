/*
 * Pointer Coordinates plug-in for Stellarium
 *
 * Copyright (C) 2014 Alexander Wolf
 * Copyright (C) 2016 Georg Zotti (Constellation code)
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "LandscapeMgr.hpp"
#include "StelMainScriptAPI.hpp"
#include "StelProjector.hpp"
#include "StelPainter.hpp"
#include "StelApp.hpp"
#include "StelCore.hpp"
#include "StelMainView.hpp"
#include "SkyGui.hpp"
#include "StelLocaleMgr.hpp"
#include "StelModuleMgr.hpp"
#include "StelFileMgr.hpp"
#include "StelGui.hpp"
#include "StelGuiItems.hpp"
#include "StelObjectMgr.hpp"
#include "StelUtils.hpp"
#include "SolarSystem.hpp"
#include "PointerCoordinates.hpp"
#include "PointerCoordinatesWindow.hpp"

#include <QFontMetrics>
#include <QSettings>
#include <QPixmap>
#include <QPair>
#include <QMetaEnum>
#include <cmath>

StelModule* PointerCoordinatesStelPluginInterface::getStelModule() const
{
	return new PointerCoordinates();
}

StelPluginInfo PointerCoordinatesStelPluginInterface::getPluginInfo() const
{
	// Allow to load the resources when used as a static plugin
	Q_INIT_RESOURCE(PointerCoordinates);

	StelPluginInfo info;
	info.id = "PointerCoordinates";
	info.displayedName = N_("Pointer Coordinates");
	info.authors = "Alexander Wolf";
	info.contact = STELLARIUM_URL;
	info.description = N_("This plugin shows the coordinates of the mouse pointer.");
	info.version = POINTERCOORDINATES_PLUGIN_VERSION;
	info.license = POINTERCOORDINATES_PLUGIN_LICENSE;
	return info;
}

PointerCoordinates::PointerCoordinates()
	: currentPlace(TopRight)
	, currentCoordinateSystem(RaDecJ2000)
	, flagShowCoordinates(false)
	, flagEnableAtStartup(false)
	, flagShowCoordinatesButton(false)
	, flagShowConstellation(false)
	, flagShowCrossedLines(false)
	, flagShowSkybright(true) // Sky luminance mod. TODO: hook up to configuration.
	, textColor(Vec3f(1,0.5,0))
	, coordinatesPoint(Vec3d(0,0,0))
	, fontSize(14)
	, toolbarButton(Q_NULLPTR)
{
	setObjectName("PointerCoordinates");
	mainWindow = new PointerCoordinatesWindow();
	StelApp &app = StelApp::getInstance();
	conf = app.getSettings();
	gui = dynamic_cast<StelGui*>(app.getGui());

}

PointerCoordinates::~PointerCoordinates()
{
	delete mainWindow;
}

void PointerCoordinates::init()
{
	if (!conf->childGroups().contains("PointerCoordinates"))
	{
		qDebug() << "[PointerCoordinates] no coordinates section exists in main config file - creating with defaults";
		restoreDefaultConfiguration();
	}

	// populate settings from main config file.
	loadConfiguration();

	addAction("actionShow_MousePointer_Coordinates",        N_("Pointer Coordinates"), N_("Show coordinates of the mouse pointer"), "enabled", "");
	addAction("actionShow_MousePointer_Coordinates_dialog", N_("Pointer Coordinates"), N_("Show settings dialog"), mainWindow, "visible");

	connect(StelApp::getInstance().getCore(), SIGNAL(configurationDataSaved()), this, SLOT(saveSettings()));

	enableCoordinates(getFlagEnableAtStartup());
	setFlagShowCoordinatesButton(flagShowCoordinatesButton);
	setFlagShowConstellation(flagShowConstellation);
	setFlagShowCrossedLines(flagShowCrossedLines);

}


void PointerCoordinates::draw(StelCore *core)
{
	if (!isEnabled())
		return;

	const StelProjectorP prj = core->getProjection(StelCore::FrameJ2000, StelCore::RefractionAuto);
	StelPainter sPainter(prj);
	sPainter.setColor(textColor, 1.f);
	font.setPixelSize(getFontSize());
	sPainter.setFont(font);

	Vec3d mousePosition = core->getMouseJ2000Pos();

	bool withDecimalDegree = StelApp::getInstance().getFlagShowDecimalDegrees();
	bool useSouthAzimuth = StelApp::getInstance().getFlagSouthAzimuthUsage();
	StelProjector::StelProjectorParams params = core->getCurrentStelProjectorParams();

	QString coordsSystem, cxt, cyt;
	double cx, cy;
	float ppx = static_cast<float>(params.devicePixelsPerPixel);
	int x, y;
	switch (getCurrentCoordinateSystem())
	{
		case RaDecJ2000:
		{
			StelUtils::rectToSphe(&cx,&cy,mousePosition); // Calculate RA/DE (J2000.0) and show it...
			coordsSystem = qc_("RA/Dec (J2000.0)", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(cx, 5, false, true);
				cyt = StelUtils::radToDecDegStr(cy);
			}
			else
			{
				cxt = StelUtils::radToHmsStr(cx, true);
				cyt = StelUtils::radToDmsStr(cy, true);
			}
			break;
		}
		case RaDec:
		{
			StelUtils::rectToSphe(&cx,&cy,core->j2000ToEquinoxEqu(mousePosition, StelCore::RefractionOff)); // Calculate RA/DE and show it...
			coordsSystem = qc_("RA/Dec", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(cx, 5, false, true);
				cyt = StelUtils::radToDecDegStr(cy);
			}
			else
			{
				cxt = StelUtils::radToHmsStr(cx, true);
				cyt = StelUtils::radToDmsStr(cy, true);
			}
			break;
		}
		case AltAzi:
		{
			StelUtils::rectToSphe(&cy,&cx,core->j2000ToAltAz(mousePosition, StelCore::RefractionAuto));
			const double direction = (useSouthAzimuth ? 2. : 3.); // N is zero, E is 90 degrees
			cy = direction*M_PI - cy;
			if (cy > M_PI*2)
				cy -= M_PI*2;

			coordsSystem = qc_("Az/Alt", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(cy);
				cyt = StelUtils::radToDecDegStr(cx);
			}
			else
			{
				cxt = StelUtils::radToDmsStr(cy);
				cyt = StelUtils::radToDmsStr(cx);
			}
			break;
		}
		case Galactic:
		{
			StelUtils::rectToSphe(&cx,&cy,core->j2000ToGalactic(mousePosition)); // Calculate galactic position and show it...
			coordsSystem = qc_("Gal. Long/Lat", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(cx);
				cyt = StelUtils::radToDecDegStr(cy);
			}
			else
			{
				cxt = StelUtils::radToDmsStr(cx, true);
				cyt = StelUtils::radToDmsStr(cy, true);
			}
			break;
		}
		case Supergalactic:
		{
			StelUtils::rectToSphe(&cx,&cy,core->j2000ToSupergalactic(mousePosition)); // Calculate supergalactic position and show it...
			coordsSystem = qc_("Supergal. Long/Lat", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(cx);
				cyt = StelUtils::radToDecDegStr(cy);
			}
			else
			{
				cxt = StelUtils::radToDmsStr(cx, true);
				cyt = StelUtils::radToDmsStr(cy, true);
			}
			break;
		}
		case Ecliptic:
		{
			double lambda, beta;
			StelUtils::rectToSphe(&cx,&cy,core->j2000ToEquinoxEqu(mousePosition, StelCore::RefractionOff));
			StelUtils::equToEcl(cx, cy, GETSTELMODULE(SolarSystem)->getEarth()->getRotObliquity(core->getJDE()), &lambda, &beta); // Calculate ecliptic position and show it...
			if (lambda<0) lambda+=2.0*M_PI;
			coordsSystem = qc_("Ecl. Long/Lat", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(lambda);
				cyt = StelUtils::radToDecDegStr(beta);
			}
			else
			{
				cxt = StelUtils::radToDmsStr(lambda, true);
				cyt = StelUtils::radToDmsStr(beta, true);
			}
			break;
		}
		case EclipticJ2000:
		{
			double lambda, beta;
			StelUtils::rectToSphe(&cx,&cy, mousePosition);
			StelUtils::equToEcl(cx, cy, GETSTELMODULE(SolarSystem)->getEarth()->getRotObliquity(2451545.0), &lambda, &beta); // Calculate ecliptic position and show it...
			if (lambda<0) lambda+=2.0*M_PI;
			coordsSystem = qc_("Ecl. Long/Lat (J2000.0)", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				cxt = StelUtils::radToDecDegStr(lambda);
				cyt = StelUtils::radToDecDegStr(beta);
			}
			else
			{
				cxt = StelUtils::radToDmsStr(lambda, true);
				cyt = StelUtils::radToDmsStr(beta, true);
			}
			break;
		}
		case HourAngle:
		{
			Vec3d v = core->j2000ToAltAz(mousePosition, StelCore::RefractionAuto);
			StelUtils::rectToSphe(&cx,&cy,Mat4d::zrotation(-core->getLocalSiderealTime())*core->altAzToEquinoxEqu(v, StelCore::RefractionOff));
			cx = 2.*M_PI-cx;
			coordsSystem = qc_("HA/Dec", "abbreviated in the plugin");
			if (withDecimalDegree)
			{
				double ha_sidereal = cx*12/M_PI;
				if (ha_sidereal>24.)
					ha_sidereal -= 24.;
				cxt = QString("%1h").arg(ha_sidereal, 0, 'f', 5);
				cyt = StelUtils::radToDecDegStr(cy);
			}
			else
			{
				cxt = StelUtils::radToHmsStr(cx);
				cyt = StelUtils::radToDmsStr(cy);
			}
			break;		
		}
	}

	QString constel;
	if (flagShowConstellation)
	{
		constel=QString(" (%1)").arg(core->getIAUConstellation(core->j2000ToEquinoxEqu(mousePosition)));
	}

	// Start sky luminance mod
	QString lumiStr;
	QString poluStr;
	if (flagShowSkybright)
	{
		// Copied from LandscapeMgr.cpp
		float latitude = core->getCurrentLocation().latitude;
		float altitude = static_cast<float>(core->getCurrentLocation().altitude);
		float temperature = 15.f;
		float relativeHumidity = 40.f;
		float lumi = 0.f;
		bool atmosphereNoScatter = false;
		float eclipseFactor = 1.f;
		//float lightPollutionLuminance = 0.f;

		SolarSystem* ssystem = static_cast<SolarSystem*>(StelApp::getInstance().getModuleMgr().getModule("SolarSystem"));
		// Compute the moon position in local coordinate
		Vec3d moonPos = ssystem->getMoon()->getAltAzPosAuto(core);
		float lunarPhaseAngle=static_cast<float>(ssystem->getMoon()->getPhaseAngle(ssystem->getEarth()->getHeliocentricEclipticPos()));
		float lunarMagnitude=ssystem->getMoon()->getVMagnitudeWithExtinction(core);

		// Compute the sun position in local coordinate
		Vec3d _sunPos = ssystem->getSun()->getAltAzPosAuto(core);

		// LP:1673283 no lunar brightening if not on Earth!
		if (core->getCurrentLocation().planetName != "Earth")
		{
			moonPos=_sunPos;
			lunarPhaseAngle=0.0f;
		}

		// Copied from Atmosphere.cpp
		if (qIsNaN(_sunPos.length()))
			_sunPos.set(0.,0.,-1.*AU);
		if (qIsNaN(moonPos.length()))
			moonPos.set(0.,0.,-1.*AU);

		// Update the eclipse intensity factor to apply on atmosphere model
		// these are for radii
		const double sun_angular_size = atan(696000./AU/_sunPos.length());
		const double moon_angular_size = atan(1738./AU/moonPos.length());
		const double touch_angle = sun_angular_size + moon_angular_size;

		// determine luminance falloff during solar eclipses
		_sunPos.normalize();
		moonPos.normalize();
		// Calculate the atmosphere RGB for each point of the grid. We can use abbreviated numbers here.
		Vec3f sunPosF=_sunPos.toVec3f();
		Vec3f moonPosF=moonPos.toVec3f();

		double separation_angle = std::acos(_sunPos.dot(moonPos));  // angle between them
		// qDebug("touch at %f\tnow at %f (%f)\n", touch_angle, separation_angle, separation_angle/touch_angle);
		// bright stars should be visible at total eclipse
		// TODO: correct for atmospheric diffusion
		// TODO: use better coverage function (non-linear)
		// because of above issues, this algorithm darkens more quickly than reality
		// Note: On Earth only, else moon would brighten other planets' atmospheres (LP:1673283)
		if ((core->getCurrentLocation().planetName=="Earth") && (separation_angle < touch_angle))
		{
			double dark_angle = moon_angular_size - sun_angular_size;
			double min = 0.0025; // 0.005f; // 0.0001f;  // so bright stars show up at total eclipse
			if (dark_angle < 0.)
			{
				// annular eclipse
				double asun = sun_angular_size*sun_angular_size;
				min = (asun - moon_angular_size*moon_angular_size)/asun;  // minimum proportion of sun uncovered
				dark_angle *= -1;
			}

			if (separation_angle < dark_angle)
				eclipseFactor = static_cast<float>(min);
			else
				eclipseFactor = static_cast<float>(min + (1.0-min)*(separation_angle-dark_angle)/(touch_angle-dark_angle));
		}
		else
			eclipseFactor = 1.f;
		// TODO: compute eclipse factor also for Lunar eclipses! (lp:#1471546)

		skyb.setLocation(latitude * M_PI_180f, altitude, temperature, relativeHumidity);
		skyb.setSunMoon(moonPosF[2], sunPosF[2]);

		// Calculate the date from the julian day.
		int year, month, day;
		StelUtils::getDateFromJulianDay(core->getJDE(), &year, &month, &day);
		skyb.setDate(year, month, lunarPhaseAngle, lunarMagnitude);

		// mousePosition is an "unProjected" Vec3d obteined in StelCore.cpp line 2729.
		// TODO: test this.
		Vec3f pointF = mousePosition.toVec3f();
		if(!atmosphereNoScatter)
		{
			// Use mirroring for sun only
			if (pointF[2]<=0.f)
			{
				pointF[2] *= -1.f;
				moonPosF[2] *= -1.f;
				// The sky below the ground is the symmetric of the one above :
				// it looks nice and gives proper values for brightness estimation
				// Use the Skybright.cpp 's models for brightness which gives better results.
			}
			// Commenting this out for now for simplicity.
			// (sky.getFlagSchaefer())
			//{
				lumi = skyb.getLuminance(moonPosF[0]*pointF[0]+moonPosF[1]*pointF[1]+moonPosF[2]*pointF[2],
							 sunPosF[0] *pointF[0]+sunPosF[1] *pointF[1]+sunPosF[2] *pointF[2],
							 pointF[2]);
			//}
			//else // Experimental: Re-allow CIE/Preetham brightness instead.
			//{
			//	skylight.pos[0]=pointF.v[0];
			//	skylight.pos[1]=pointF.v[1];
			//	skylight.pos[2]=pointF.v[2];
			//	sky.getxyYValuev(skylight);
			//	lumi=skylight.color[2];
			//}
		}

		lumi *= eclipseFactor;
		// Add star background luminance
		lumi += 0.0001f;
		// Multiply by the input scale of the ToneConverter (is not done automatically by the xyYtoRGB method called later)
		//lumi*=eye->getInputScale();

		// TODO:
		// Perhaps better to use the Atmosphere class directly rather than copy and paste.
		// Add the light pollution luminance AFTER the scaling to avoid scaling it because it is the cause
		// of the scaling itself

		// Bortle scale is managed by SkyDrawer
		int bortleIndex = StelMainScriptAPI::getBortleScaleIndex();

		// Returns some representative value of zenith luminance in cd/m² for the given Bortle scale index.
		auto lightPollutionLuminance = StelCore::bortleScaleIndexToLuminance(bortleIndex);
		auto nelm = StelCore::luminanceToNELM(lightPollutionLuminance);
		nelm = std::round(nelm*10)*0.1; // copied from LightPollutionWidget

		// As far as I understand, fader is a 0..1 crossfade visual effect and it is safe to assume 1.0 when disabled and showing the sky.
		//lumi += fader.getInterstate()*lightPollutionLuminance;

		lumiStr = QString::number(lumi, 'e');
		poluStr = QString("Pollution lum.: %1 cd/m2, Bortle class: %2, NELM: %3").arg(QString::number(lightPollutionLuminance, 'e')).arg(bortleIndex).arg(nelm);

	}
	// End sky luminance mod


	QString coordsText = QString("%1: %2/%3%4").arg(coordsSystem, cxt, cyt, constel);
	x = getCoordinatesPlace(coordsText).first;
	y = getCoordinatesPlace(coordsText).second;
	if (getCurrentCoordinatesPlace()!=Custom)
	{
		x *= ppx;
		y *= ppx;
	}
	sPainter.drawText(x, y, coordsText);

	// Start sky luminance mod
	QFontMetrics fm(font);
	QSize fs = fm.size(Qt::TextSingleLine, coordsText);
	sPainter.drawText(x, y-fs.height(), QString("Luminance: %1 lux.").arg(lumiStr));
	sPainter.drawText(x, y-2*fs.height(), poluStr);
	// End sky luminance mod


	if (flagShowCrossedLines)
	{
		QPoint m = StelMainView::getInstance().getMousePos();
		sPainter.drawLine2d(m.x()*ppx, 0, m.x()*ppx, params.viewportXywh[3]*ppx);
		sPainter.drawLine2d(0, (params.viewportXywh[3]-m.y())*ppx, params.viewportXywh[2]*ppx, (params.viewportXywh[3]-m.y())*ppx);
	}

}

void PointerCoordinates::enableCoordinates(bool b)
{
	if (b!=flagShowCoordinates)
	{
		flagShowCoordinates = b;
		emit flagCoordinatesVisibilityChanged(b);
	}
}

double PointerCoordinates::getCallOrder(StelModuleActionName actionName) const
{
	if (actionName==StelModule::ActionDraw)
		return StelApp::getInstance().getModuleMgr().getModule("LabelMgr")->getCallOrder(actionName)+110.;
	return 0;
}

bool PointerCoordinates::configureGui(bool show)
{
	if (show)
	{
		mainWindow->setVisible(true);
	}

	return true;
}

void PointerCoordinates::restoreDefaultConfiguration(void)
{
	// Remove the whole section from the configuration file
	conf->remove("PointerCoordinates");
	// Load the default values...
	loadConfiguration();
	// ... then save them.
	saveConfiguration();
	// But this doesn't save the color, so...
	conf->beginGroup("PointerCoordinates");
	conf->setValue("text_color", "1,0.5,0");
	conf->endGroup();
}

void PointerCoordinates::loadConfiguration(void)
{
	conf->beginGroup("PointerCoordinates");

	setFlagEnableAtStartup(conf->value("enable_at_startup", false).toBool());
	textColor = Vec3f(conf->value("text_color", "1,0.5,0").toString());
	setFontSize(conf->value("font_size", 14).toInt());
	flagShowCoordinatesButton = conf->value("flag_show_button", true).toBool();	
	setCurrentCoordinatesPlaceKey(conf->value("current_displaying_place", "TopRight").toString());
	setCurrentCoordinateSystemKey(conf->value("current_coordinate_system", "RaDecJ2000").toString());
	QStringList cc = conf->value("custom_coordinates", "1,1").toString().split(",");
	setCustomCoordinatesPlace(cc[0].toInt(), cc[1].toInt());
	flagShowConstellation = conf->value("flag_show_constellation", false).toBool();
	flagShowCrossedLines = conf->value("flag_show_crossed_lines", false).toBool();

	conf->endGroup();
}

void PointerCoordinates::saveConfiguration(void)
{
	conf->beginGroup("PointerCoordinates");

	conf->setValue("enable_at_startup", getFlagEnableAtStartup());
	conf->setValue("flag_show_button", getFlagShowCoordinatesButton());
	conf->setValue("current_displaying_place", getCurrentCoordinatesPlaceKey());
	conf->setValue("current_coordinate_system", getCurrentCoordinateSystemKey());
	QPair<int, int> cc = getCustomCoordinatesPlace();
	conf->setValue("custom_coordinates", QString("%1,%2").arg(cc.first).arg(cc.second));
	//conf->setValue("text_color", "1,0.5,0");
	conf->setValue("font_size", getFontSize());
	conf->setValue("flag_show_constellation", getFlagShowConstellation());
	conf->setValue("flag_show_crossed_lines", getFlagShowCrossedLines());

	conf->endGroup();
}

void PointerCoordinates::setFlagShowCoordinatesButton(bool b)
{
	if (gui!=Q_NULLPTR)
	{
		if (b==true) {
			if (toolbarButton==Q_NULLPTR) {
				// Create the button
				toolbarButton = new StelButton(Q_NULLPTR,
							       QPixmap(":/PointerCoordinates/bt_PointerCoordinates_On.png"),
							       QPixmap(":/PointerCoordinates/bt_PointerCoordinates_Off.png"),
							       QPixmap(":/graphicGui/miscGlow32x32.png"),
							       "actionShow_MousePointer_Coordinates",
							       false,
							       "actionShow_MousePointer_Coordinates_dialog");
			}
			gui->getButtonBar()->addButton(toolbarButton, "065-pluginsGroup");
		} else {
			gui->getButtonBar()->hideButton("actionShow_MousePointer_Coordinates");
		}
	}
	flagShowCoordinatesButton = b;
}

// Set the current place of the string with coordinates from its key
void PointerCoordinates::setCurrentCoordinatesPlaceKey(QString key)
{
	const QMetaEnum& en = metaObject()->enumerator(metaObject()->indexOfEnumerator("CoordinatesPlace"));
	CoordinatesPlace coordPlace = static_cast<CoordinatesPlace>(qMax(0, en.keyToValue(key.toLatin1().data())));
	setCurrentCoordinatesPlace(coordPlace);
}

//Get the current place of the string with coordinates
QString PointerCoordinates::getCurrentCoordinatesPlaceKey() const
{
	return metaObject()->enumerator(metaObject()->indexOfEnumerator("CoordinatesPlace")).key(currentPlace);
}

void PointerCoordinates::setCurrentCoordinateSystemKey(QString key)
{
	const QMetaEnum& en = metaObject()->enumerator(metaObject()->indexOfEnumerator("CoordinateSystem"));
	CoordinateSystem coordSystem = static_cast<CoordinateSystem>(qMax(0, en.keyToValue(key.toLatin1().data())));
	setCurrentCoordinateSystem(coordSystem);
}

QString PointerCoordinates::getCurrentCoordinateSystemKey() const
{
	return metaObject()->enumerator(metaObject()->indexOfEnumerator("CoordinateSystem")).key(currentCoordinateSystem);
}

QPair<int, int> PointerCoordinates::getCoordinatesPlace(QString text)
{
	int x = 0, y = 0;
	static const float coeff = 1.5;
	QFontMetrics fm(font);
	QSize fs = fm.size(Qt::TextSingleLine, text);
	switch(getCurrentCoordinatesPlace())
	{
		case TopCenter:
		{
			x = gui->getSkyGui()->getSkyGuiWidth()/2 - fs.width()/2;
			y = gui->getSkyGui()->getSkyGuiHeight() - static_cast<int>(fs.height()*coeff);
			break;
		}
		case TopRight:
		{
			x = 3*gui->getSkyGui()->getSkyGuiWidth()/4 - fs.width()/2;
			y = gui->getSkyGui()->getSkyGuiHeight() - static_cast<int>(fs.height()*coeff);
			break;
		}
		case RightBottomCorner:
		{
			x = gui->getSkyGui()->getSkyGuiWidth() - static_cast<int>(fs.width() + 10*coeff);
			y = fs.height();
			break;
		}
		case NearMouseCursor:
		{
			QPoint m = StelMainView::getInstance().getMousePos();
			x = m.x() + 3;
			y = gui->getSkyGui()->getSkyGuiHeight() - m.y() + 5;
			break;
		}
		case Custom:
		{
			QPair<int, int> xy = getCustomCoordinatesPlace();
			x = xy.first;
			y = gui->getSkyGui()->getSkyGuiHeight() - xy.second - fs.height();
			break;
		}
	}
	return qMakePair(x, y);
}

void PointerCoordinates::setCustomCoordinatesPlace(int x, int y)
{
	customPosition = qMakePair(x, y);
}


