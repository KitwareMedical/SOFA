/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "SofaVideoRecorderManager.h"

#include <iostream>
#ifndef SOFA_QT4
#include <qlineedit.h>
#include <qcombobox.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qlayout.h>
#include <qradiobutton.h>
#endif

namespace sofa
{
namespace gui
{
namespace qt
{

CaptureOptionsWidget::CaptureOptionsWidget( QWidget * parent)
    : QWidget(parent)
{

    QVBoxLayout *layout=new QVBoxLayout(this);

    QHBoxLayout *HLayoutFramerate = new QHBoxLayout();
    QLabel *labelFramerate=new QLabel(QString("Framerate (in img/s): "), this);
    framerateSpinBox = new QSpinBox(this);
    framerateSpinBox->setMinValue(1);
    framerateSpinBox->setMaxValue(120);
    framerateSpinBox->setValue(60);
    HLayoutFramerate->addWidget (labelFramerate);
    HLayoutFramerate->addWidget (framerateSpinBox);

    realtimeCheckBox = new QCheckBox(QString("Real-Time recording"), this);

    QHBoxLayout *HLayoutFrameskip = new QHBoxLayout();
    QLabel *labelFrameskip=new QLabel(QString("Skip frames before capture (fast replay): "), this);
    frameskipSpinBox = new QSpinBox(this);
    frameskipSpinBox->setMinValue(0);
    frameskipSpinBox->setMaxValue(100);
    frameskipSpinBox->setValue(0);
    HLayoutFrameskip->addWidget (labelFrameskip);
    HLayoutFrameskip->addWidget (frameskipSpinBox);

    layout->addLayout(HLayoutFramerate);
    layout->addWidget(realtimeCheckBox);
    layout->addLayout(HLayoutFrameskip);

    //this->addLayout(layout);
}

MovieOptionsWidget::MovieOptionsWidget( QWidget * parent)
    : QWidget(parent)
{
    //Build codec list
    listCodecs.push_back(Codec("mpeg", "Mpeg1 (Bad quality but readable everywhere)"));
    listCodecs.push_back(Codec("mp4", "MP4/Mpeg4 (Good ratio visual quality/bitrate, good compatibility)"));
    listCodecs.push_back(Codec("mp4","h264", "MP4/H264 (Best ratio visual quality/bitrate, requires libx264)"));
    listCodecs.push_back(Codec("avi","lossless", "Lossless (No loss of information, best for post-processing and re-encodings)"));

    QVBoxLayout *layout=new QVBoxLayout(this);

    QHBoxLayout *HLayoutCodec = new QHBoxLayout();
    QLabel *labelCodec=new QLabel(QString("Codec: "), this);
    codecComboBox = new QComboBox(this);
    for(unsigned int i=0; i<listCodecs.size(); i++)
        codecComboBox->insertItem(QString(listCodecs[i].description.c_str()));
    codecComboBox->setCurrentIndex(2);
    HLayoutCodec->addWidget (labelCodec);
    HLayoutCodec->addWidget (codecComboBox);

    QHBoxLayout *HLayoutBitrate = new QHBoxLayout();
    QLabel *labelBitrate=new QLabel(QString("Bitrate (in KB/s): "), this);
    bitrateSpinBox = new QSpinBox(this);
    bitrateSpinBox->setMinValue(100);
    bitrateSpinBox->setMaxValue(40960);
    bitrateSpinBox->setValue(5000);
    HLayoutBitrate->addWidget (labelBitrate);
    HLayoutBitrate->addWidget (bitrateSpinBox);

    layout->addLayout(HLayoutCodec);
    layout->addLayout(HLayoutBitrate);

    //this->addLayout(layout);
}

SofaVideoRecorderManager::SofaVideoRecorderManager()
{
    setupUi(this);
    captureOptionsWidget = new CaptureOptionsWidget(this);
    movieOptionsWidget = new MovieOptionsWidget(this);

    internalAddWidget(VideoRecorderOptionGroupBox, captureOptionsWidget);
    internalAddWidget(VideoRecorderOptionGroupBox, movieOptionsWidget);

#ifdef SOFA_HAVE_FFMPEG
    MovieRecordingTypeRadioButton->setChecked(true);
#else
    MovieRecordingTypeRadioButton->setHidden(true);
#endif
    onChangeRecordingType();
}


std::string SofaVideoRecorderManager::getCodecExtension()
{
    unsigned int index = movieOptionsWidget->codecComboBox->currentItem();
    return movieOptionsWidget->listCodecs[index].extension;
}

std::string SofaVideoRecorderManager::getCodecName()
{
    unsigned int index = movieOptionsWidget->codecComboBox->currentItem();
    return movieOptionsWidget->listCodecs[index].codec;
}

unsigned int SofaVideoRecorderManager::getFramerate()
{
    return captureOptionsWidget->framerateSpinBox->value();
}

unsigned int SofaVideoRecorderManager::getBitrate()
{
    return movieOptionsWidget->bitrateSpinBox->value()*1024;
}

bool SofaVideoRecorderManager::realtime()
{
    return captureOptionsWidget->realtimeCheckBox->isChecked();
}

unsigned int SofaVideoRecorderManager::getFrameskip()
{
    return captureOptionsWidget->frameskipSpinBox->value();
}


void SofaVideoRecorderManager::updateContent()
{
    movieOptionsWidget->setHidden(currentRecordingType != MOVIE);
}

void SofaVideoRecorderManager::onChangeRecordingType()
{
    currentRecordingType = (MovieRecordingTypeRadioButton->isChecked()) ? MOVIE : SCREENSHOTS;

    updateContent();
}

void SofaVideoRecorderManager::internalAddWidget(QWidget* parent, QWidget* widgetToAdd)
{

#ifdef SOFA_QT4
    parent->layout()->addWidget(widgetToAdd);
#else
    parent->layout()->add(widgetToAdd);
#endif
}

void SofaVideoRecorderManager::close()
{
    this->hide();
}

} //namespace qt

} //namespace gui

}//namespace sofa
