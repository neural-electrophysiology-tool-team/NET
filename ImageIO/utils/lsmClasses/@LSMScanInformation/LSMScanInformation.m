classdef LSMScanInformation
  %LSMSCANINFORMATION Class representation of ScanInformation in LSM files
  % The field u32OffsetScanInformation of the "CZ-Private tag" contains the
  % file offset to the start of a block with information of the device settings
  % used during acquisition. Note that the image size, channel number and pixel
  % distance of the image stored in the file can be different form the settings
  % used during acquisition because the image could have been modified by offline
  % operations. Information about the image size, channel number and pixel distance
  % of the image contents stored in the file can be found in the fields of the "CZ-Private tag"
  %
  % AUTHOR: Stefano Masneri
  % Date: 14.3.2017
  
  
  properties
    entries;        % Cell containing all the properties read from the ScanInformation

    propertyList;   % Cell array matching hex codes to metadata, based on LSM 
                    % file format specifications
  end
  
  properties (Hidden = true)
    idx = 1;
  end
  
  methods
    function obj = LSMScanInformation(lsmPtr, byteOrder)
    %LSMOVERLAY Constructor
    % Assumes the file pointer in the correct position already
      obj = obj.createPropertyList();
      
      obj.entries = cell(207, 2);
      
      % The algorithm for reading the scaninfo database depends on keeping track of a
      % "level" hierarchy.  As the database is read, some entries are level
      % instructions, navigating up or down the hierarchy.  The database is done when
      % the hierarchy steps back to level 0
      level = 0;
      while (1)
        
        taghex      = dec2hex(fread(lsmPtr, 1, 'uint32', byteOrder));
        typecode    = fread(lsmPtr, 1, 'uint32', byteOrder);
        size        = fread(lsmPtr, 1, 'uint32', byteOrder);
        tag         = obj.getTag(['h' taghex]);
        
        switch typecode
          case 0 % Special case: this is a level instruction entry
            if (hex2dec(taghex) == hex2dec('FFFFFFFF'))
              level = level - 1;
            else
              level = level + 1;
            end
          case 2 % string
            count = size;
            value = char(fread(lsmPtr, count, 'uchar', byteOrder)');
            value = value(1:end-1);
          case 4 % int32
            count = size / 4;
            value = fread(lsmPtr, count, 'uint32', byteOrder);
          case 5 % float64
            count = size / 8;
            value = fread(lsmPtr, count, 'float64', byteOrder);
          otherwise
            error('LSMScanInformation: error parsing LSM file')
        end
        
        % If this was just a level instruction entry ignore it, otherwise try to
        % record entry.
        if typecode > 0
          if isempty(tag)
            propName = [ 'h' taghex ];
          else
            propName = tag;
          end
          obj = obj.addEntry(propName, value);
        end
        
        if 0 == level
          break;
        end
      end
    end
    
    function [keys, values] = getAllEntries(obj)
      keys = obj.entries.keys();
      values = obj.entries.values();
    end
  end
  
  methods (Access = private)
    
    function tag = getTag(obj, tagHex)
      tag = [];
      for k = 1:length(obj.propertyList)
        if strcmp(tagHex, obj.propertyList{k}{1})
          tag = obj.propertyList{k}{2};
          break;
        end
      end
    end
    
    function obj = addEntry(obj, propName, value)
      try
        obj.entries{obj.idx, 1} = propName;
        obj.entries{obj.idx, 2} = value;
        obj.idx = obj.idx + 1;
      catch
        error('LSMScanInformation: error adding field to entries')
      end
    end
    
    function obj = createPropertyList(obj)
      cnt = 1;
      scan{cnt}  = {'h10000000',    'RECORDINGS'                      }; cnt = cnt+1;
      scan{cnt}  = {'h10000001',    'ENTRY_NAME'                      }; cnt = cnt+1;
      scan{cnt}  = {'h10000002',    'ENTRY_DESCRIPTION'               }; cnt = cnt+1;
      scan{cnt}  = {'h10000003',    'ENTRY_NOTES'                   	}; cnt = cnt+1;
      scan{cnt}  = {'h10000004',    'ENTRY_OBJECTIVE'                 }; cnt = cnt+1;
      scan{cnt}  = {'h10000005',    'PROCESSING_SUMMARY'              }; cnt = cnt+1;
      scan{cnt}  = {'h10000006',    'SPECIAL_SCAN'                    }; cnt = cnt+1;
      scan{cnt}  = {'h10000007',    'SCAN_TYPE'                       }; cnt = cnt+1;
      scan{cnt}  = {'h10000008',    'SCAN_MODE'                       }; cnt = cnt+1;
      scan{cnt}  = {'h10000009',    'STACKS_COUNT'                    }; cnt = cnt+1;
      scan{cnt}  = {'h1000000A',    'LINES_PER_PLANE'                 }; cnt = cnt+1;
      scan{cnt}  = {'h1000000B',    'SAMPLES_PER_LINE'                }; cnt = cnt+1;
      scan{cnt}  = {'h1000000C',    'PLANES_PER_VOLUME'               }; cnt = cnt+1;
      scan{cnt}  = {'h1000000D',    'IMAGES_WIDTH'                    }; cnt = cnt+1;
      scan{cnt}  = {'h1000000E',    'IMAGES_HEIGHT'                   }; cnt = cnt+1;
      scan{cnt}  = {'h1000000F',    'NUMBER_OF_PLANES'                }; cnt = cnt+1;
      scan{cnt}  = {'h10000010',    'IMAGES_NUMBER_STACKS'            }; cnt = cnt+1;
      scan{cnt}  = {'h10000011',    'IMAGES_NUMBER_CHANNELS'          }; cnt = cnt+1;
      scan{cnt}  = {'h10000012',    'LINESCAN_X_Y'                    }; cnt = cnt+1;
      scan{cnt}  = {'h10000013',    'SCAN_DIRECTION'                  }; cnt = cnt+1;
      scan{cnt}  = {'h10000014',    'TIME_SERIES'                   	}; cnt = cnt+1;
      scan{cnt}  = {'h10000015',    'ORIGINAL_SCAN_DATA'             	}; cnt = cnt+1;
      scan{cnt}  = {'h10000016',    'ZOOM_X'                          }; cnt = cnt+1;
      scan{cnt}  = {'h10000017',    'ZOOM_Y'                          }; cnt = cnt+1;
      scan{cnt}  = {'h10000018',    'ZOOM_Z'                          }; cnt = cnt+1;
      scan{cnt}  = {'h10000019',    'SAMPLE_0_X'                      }; cnt = cnt+1;
      scan{cnt}  = {'h1000001A',    'SAMPLE_0_Y'                      }; cnt = cnt+1;
      scan{cnt}  = {'h1000001B',    'SAMPLE_0_Z'                      }; cnt = cnt+1;
      scan{cnt}  = {'h1000001C',    'SAMPLE_SPACING'                  }; cnt = cnt+1;
      scan{cnt}  = {'h1000001D',    'LINE_SPACING'                    }; cnt = cnt+1;
      scan{cnt}  = {'h1000001E',    'PLANE_SPACING'                   }; cnt = cnt+1;
      scan{cnt}  = {'h1000001F',    'PLANE_WIDTH'                     }; cnt = cnt+1;
      scan{cnt}  = {'h10000020',    'PLANE_HEIGHT'                    }; cnt = cnt+1;
      scan{cnt}  = {'h10000021',    'VOLUME_DEPTH'                  	}; cnt = cnt+1;
      scan{cnt}  = {'h10000034',    'ROTATION'                      	}; cnt = cnt+1;
      scan{cnt}  = {'h10000035',    'PRECESSION'                    	}; cnt = cnt+1;
      scan{cnt}  = {'h10000036',    'SAMPLE_0_TIME'                   }; cnt = cnt+1;
      scan{cnt}  = {'h10000037',    'START_SCAN_TRIGGER_IN'         	}; cnt = cnt+1;
      scan{cnt}  = {'h10000038',    'START_SCAN_TRIGGER_OUT'          }; cnt = cnt+1;
      scan{cnt}  = {'h10000039',    'START_SCAN_EVENT'                }; cnt = cnt+1;
      scan{cnt}  = {'h10000040',    'START_SCAN_TIME'                 }; cnt = cnt+1;
      scan{cnt}  = {'h10000041',    'STOP_SCAN_TRIGGER_IN'            }; cnt = cnt+1;
      scan{cnt}  = {'h10000042',    'STOP_SCAN_TRIGGER_OUT'           }; cnt = cnt+1;
      scan{cnt}  = {'h10000043',    'STOP_SCAN_EVENT'                 }; cnt = cnt+1;
      scan{cnt}  = {'h10000044',    'START_SCAN_TIME2'                }; cnt = cnt+1;
      scan{cnt}  = {'h10000045',    'USE_ROIS'                        }; cnt = cnt+1;
      scan{cnt}  = {'h10000046',    'USE_REDUCED_MEMORY_ROIS'         }; cnt = cnt+1;
      scan{cnt}  = {'h10000047',    'USER'                            }; cnt = cnt+1;
      scan{cnt}  = {'h10000048',    'USE_BCCORECCTION'                }; cnt = cnt+1;
      scan{cnt}  = {'h10000049',    'POSITION_BCCORRECTION1'        	}; cnt = cnt+1;
      scan{cnt}  = {'h10000050',    'POSITION_BCCORRECTION2'        	}; cnt = cnt+1;
      scan{cnt}  = {'h10000051',    'INTERPOLATIONY'                  }; cnt = cnt+1;
      scan{cnt}  = {'h10000052',    'CAMERA_BINNING'                  }; cnt = cnt+1;
      scan{cnt}  = {'h10000053',    'CAMERA_SUPERSAMPLING'            }; cnt = cnt+1;
      scan{cnt}  = {'h10000054',    'CAMERA_FRAME_WIDTH'              }; cnt = cnt+1;
      scan{cnt}  = {'h10000055',    'CAMERA_FRAME_HEIGHT'             }; cnt = cnt+1;
      scan{cnt}  = {'h10000056',    'CAMERA_OFFSETX'                  }; cnt = cnt+1;
      scan{cnt}  = {'h10000057',    'CAMERA_OFFSETY'                  }; cnt = cnt+1;
      scan{cnt}  = {'h10000059',    'RT_BINNNING'                     }; cnt = cnt+1;
      scan{cnt}  = {'h1000005a',    'RT_FRAME_WIDTH'                  }; cnt = cnt+1;
      scan{cnt}  = {'h1000005b',    'RT_FRAME_HEIGHT'                 }; cnt = cnt+1;
      scan{cnt}  = {'h1000005c',    'RT_REGION_WIDTH'                 }; cnt = cnt+1;
      scan{cnt}  = {'h1000005d',    'RT_REGION_HEIGHT'                }; cnt = cnt+1;
      scan{cnt}  = {'h1000005e',    'RT_OFFSETX'                      }; cnt = cnt+1;
      scan{cnt}  = {'h1000005f',    'RT_OFFSETY'                      }; cnt = cnt+1;
      scan{cnt}  = {'h10000060',    'RT_ZOOM'                         }; cnt = cnt+1;
      scan{cnt}  = {'h10000061',    'RT_LINEPERIOD'                   }; cnt = cnt+1;
      scan{cnt}  = {'h10000062',    'PRESCAN'                         }; cnt = cnt+1;
      scan{cnt}  = {'h10000063',    'SCAN_DIRECTIONZ'                 }; cnt = cnt+1;
      scan{cnt}  = {'h11000000',    'TIMERS'                          }; cnt = cnt+1;
      scan{cnt}  = {'h12000000',    'TIMER'                           }; cnt = cnt+1;
      scan{cnt}  = {'h12000001',    'TIMER_NAME'                      }; cnt = cnt+1;
      scan{cnt}  = {'h12000003',    'INTERVAL'                        }; cnt = cnt+1;
      scan{cnt}  = {'h12000004',    'TRIGGER_IN'                      }; cnt = cnt+1;
      scan{cnt}  = {'h12000005',    'TRIGGER_OUT'                     }; cnt = cnt+1;
      scan{cnt}  = {'h13000000',    'MARKERS'                         }; cnt = cnt+1;
      scan{cnt}  = {'h14000000',    'MARKER'                          }; cnt = cnt+1;
      scan{cnt}  = {'h14000001',    'MARKER_NAME'                     }; cnt = cnt+1;
      scan{cnt}  = {'h14000002',    'DESCRIPTION'                     }; cnt = cnt+1;
      scan{cnt}  = {'h14000003',    'TRIGGER_IN'                      }; cnt = cnt+1;
      scan{cnt}  = {'h14000004',    'TRIGGER_OUT'                   	}; cnt = cnt+1;
      scan{cnt}  = {'h20000000',    'TRACKS'                          }; cnt = cnt+1;
      scan{cnt}  = {'h30000000',    'LASERS'                          }; cnt = cnt+1;
      scan{cnt}  = {'h40000000',    'TRACK'                           }; cnt = cnt+1;
      scan{cnt}  = {'h40000001',    'MULTIPLEX_TYPE'                  }; cnt = cnt+1;
      scan{cnt}  = {'h40000002',    'MULTIPLEX_ORDER'	                }; cnt = cnt+1;
      scan{cnt}  = {'h40000003',    'SAMPLING_MODE'                   }; cnt = cnt+1;
      scan{cnt}  = {'h40000004',    'SAMPLING_METHOD'	                }; cnt = cnt+1;
      scan{cnt}  = {'h40000005',    'SAMPLING_NUMBER'	                }; cnt = cnt+1;
      scan{cnt}  = {'h40000006',    'ENTRY_ACQUIRE'                   }; cnt = cnt+1;
      scan{cnt}  = {'h40000007',    'OBSERVATION_TIME'	              }; cnt = cnt+1;
      scan{cnt}  = {'h4000000B',    'TIME_BETWEEN_STACKS'	            }; cnt = cnt+1;
      scan{cnt}  = {'h4000000C',    'TRACK_NAME'                      }; cnt = cnt+1;
      scan{cnt}  = {'h4000000D',    'COLLIMATOR1_NAME'	              }; cnt = cnt+1;
      scan{cnt}  = {'h4000000E',    'COLLIMATOR1_POSITION'            }; cnt = cnt+1;
      scan{cnt}  = {'h4000000F',    'COLLIMATOR2_NAME'                }; cnt = cnt+1;
      scan{cnt}  = {'h40000010',    'COLLIMATOR2_POSITION'            }; cnt = cnt+1;
      scan{cnt}  = {'h40000011',    'BLEACH_TRACK'                    }; cnt = cnt+1;
      scan{cnt}  = {'h40000012',    'BLEACH_AFTER_SCAN_NUMBER'        }; cnt = cnt+1;
      scan{cnt}  = {'h40000013',    'BLEACH_SCAN_NUMBER'              }; cnt = cnt+1;
      scan{cnt}  = {'h40000014',    'TRIGGER_IN'                      }; cnt = cnt+1;
      scan{cnt}  = {'h40000015',    'TRIGGER_OUT'                     }; cnt = cnt+1;
      scan{cnt}  = {'h40000016',    'IS_RATIO_TRACK'                  }; cnt = cnt+1;
      scan{cnt}  = {'h40000017',    'BLEACH_COUNT'                    }; cnt = cnt+1;
      scan{cnt}  = {'h40000018',    'SPI_CENTER_WAVELENGTH'           }; cnt = cnt+1;
      scan{cnt}  = {'h40000019',    'PIXEL_TIME'                      }; cnt = cnt+1;
      scan{cnt}  = {'h40000020',    'ID_CONDENSOR_FRONTLENS'          }; cnt = cnt+1;
      scan{cnt}  = {'h40000021',    'CONDENSOR_FRONTLENS'             }; cnt = cnt+1;
      scan{cnt}  = {'h40000022',    'ID_FIELD_STOP'                   }; cnt = cnt+1;
      scan{cnt}  = {'h40000023',    'FIELD_STOP_VALUE'                }; cnt = cnt+1;
      scan{cnt} = {'h40000024',    'ID_CONDENSOR_APERTURE'            }; cnt = cnt+1;
      scan{cnt} = {'h40000025',    'CONDENSOR_APERTURE'               }; cnt = cnt+1;
      scan{cnt} = {'h40000026',    'ID_CONDENSOR_REVOLVER'            }; cnt = cnt+1;
      scan{cnt} = {'h40000027',    'CONDENSOR_FILTER'                 }; cnt = cnt+1;
      scan{cnt} = {'h40000028',    'ID_TRANSMISSION_FILTER1'          }; cnt = cnt+1;
      scan{cnt} = {'h40000029',    'ID_TRANSMISSION1'                 }; cnt = cnt+1;
      scan{cnt} = {'h40000030',    'ID_TRANSMISSION_FILTER2'          }; cnt = cnt+1;
      scan{cnt} = {'h40000031',    'ID_TRANSMISSION2'                 }; cnt = cnt+1;
      scan{cnt} = {'h40000032',    'REPEAT_BLEACH'                    }; cnt = cnt+1;
      scan{cnt} = {'h40000033',    'ENABLE_SPOT_BLEACH_POS'           }; cnt = cnt+1;
      scan{cnt} = {'h40000034',    'SPOT_BLEACH_POSX'	                }; cnt = cnt+1;
      scan{cnt} = {'h40000035',    'SPOT_BLEACH_POSY'	                }; cnt = cnt+1;
      scan{cnt} = {'h40000036',    'BLEACH_POSITION_Z'	              }; cnt = cnt+1;
      scan{cnt} = {'h40000037',    'ID_TUBELENS'                      }; cnt = cnt+1;
      scan{cnt} = {'h40000038',    'ID_TUBELENS_POSITION'          	  }; cnt = cnt+1;
      scan{cnt} = {'h40000039',    'TRANSMITTED_LIGHT'             	  }; cnt = cnt+1;
      scan{cnt} = {'h4000003a',    'REFLECTED_LIGHT'	                }; cnt = cnt+1;
      scan{cnt} = {'h4000003b',    'TRACK_SIMULTAN_GRAB_AND_BLEACH'   }; cnt = cnt+1;
      scan{cnt} = {'h4000003c',    'BLEACH_PIXEL_TIME'	              }; cnt = cnt+1;
      scan{cnt} = {'h50000000',    'LASER'	                          }; cnt = cnt+1;
      scan{cnt} = {'h50000001',    'LASER_NAME'                       }; cnt = cnt+1;
      scan{cnt} = {'h50000002',    'LASER_ACQUIRE'                    }; cnt = cnt+1;
      scan{cnt} = {'h50000003',    'LASER_POWER'                      }; cnt = cnt+1;
      scan{cnt} = {'h60000000',    'DETECTION_CHANNELS'               }; cnt = cnt+1;
      scan{cnt} = {'h70000000',    'DETECTION_CHANNEL'	              }; cnt = cnt+1;
      scan{cnt} = {'h70000003',    'DETECTOR_GAIN'                    }; cnt = cnt+1;
      scan{cnt} = {'h70000005',    'AMPLIFIER_GAIN'                   }; cnt = cnt+1;
      scan{cnt} = {'h70000007',    'AMPLIFIER_OFFSET'                 }; cnt = cnt+1;
      scan{cnt} = {'h70000009',    'PINHOLE_DIAMETER'                 }; cnt = cnt+1;
      scan{cnt} = {'h7000000B',    'ENTRY_ACQUIRE'                    }; cnt = cnt+1;
      scan{cnt} = {'h7000000C',    'DETECTOR_NAME'                    }; cnt = cnt+1;
      scan{cnt} = {'h7000000D',    'AMPLIFIER_NAME'                   }; cnt = cnt+1;
      scan{cnt} = {'h7000000E',    'PINHOLE_NAME'                     }; cnt = cnt+1;
      scan{cnt} = {'h7000000F',    'FILTER_SET_NAME'                  }; cnt = cnt+1;
      scan{cnt} = {'h70000010',    'FILTER_NAME'                   	  }; cnt = cnt+1;
      scan{cnt} = {'h70000013',    'INTEGRATOR_NAME'                  }; cnt = cnt+1;
      scan{cnt} = {'h70000014',    'DETECTION_CHANNEL_NAME'           }; cnt = cnt+1;
      scan{cnt} = {'h70000015',    'DETECTOR_GAIN_BC1'                }; cnt = cnt+1;
      scan{cnt} = {'h70000016',    'DETECTOR_GAIN_BC2'                }; cnt = cnt+1;
      scan{cnt} = {'h70000017',    'AMPLIFIER_GAIN_BC1'               }; cnt = cnt+1;
      scan{cnt} = {'h70000018',    'AMPLIFIER_GAIN_BC2'               }; cnt = cnt+1;
      scan{cnt} = {'h70000019',    'AMPLIFIER_OFFSET_BC1'             }; cnt = cnt+1;
      scan{cnt} = {'h70000020',    'AMPLIFIER_OFFSET_BC2'             }; cnt = cnt+1;
      scan{cnt} = {'h70000021',    'SPECTRAL_SCAN_CHANNELS'           }; cnt = cnt+1;
      scan{cnt} = {'h70000022',    'SPI_WAVE_LENGTH_START'         	  }; cnt = cnt+1;
      scan{cnt} = {'h70000023',    'SPI_WAVELENGTH_END'               }; cnt = cnt+1;
      scan{cnt} = {'h70000026',    'DYE_NAME'                         }; cnt = cnt+1;
      scan{cnt} = {'h70000027',    'DYE_FOLDER'                       }; cnt = cnt+1;
      scan{cnt} = {'h80000000',    'ILLUMINATION_CHANNELS'            }; cnt = cnt+1;
      scan{cnt} = {'h90000000',    'ILLUMINATION_CHANNEL'             }; cnt = cnt+1;
      scan{cnt} = {'h90000001',    'ILL_NAME'                         }; cnt = cnt+1;
      scan{cnt} = {'h90000002',    'POWER'                            }; cnt = cnt+1;
      scan{cnt} = {'h90000003',    'WAVELENGTH'                       }; cnt = cnt+1;
      scan{cnt} = {'h90000004',    'ACQUIRE'                          }; cnt = cnt+1;
      scan{cnt} = {'h90000005',    'DETCHANNEL_NAME'                  }; cnt = cnt+1;
      scan{cnt} = {'h90000006',    'POWER_BC1'                        }; cnt = cnt+1;
      scan{cnt} = {'h90000007',    'POWER_BC2'                        }; cnt = cnt+1;
      scan{cnt} = {'hA0000000',    'BEAM_SPLITTERS'                   }; cnt = cnt+1;
      scan{cnt} = {'hB0000000',    'BEAM_SPLITTER'                    }; cnt = cnt+1;
      scan{cnt} = {'hB0000001',    'FILTER_SET'                       }; cnt = cnt+1;
      scan{cnt} = {'hB0000002',    'FILTER'                           }; cnt = cnt+1;
      scan{cnt} = {'hB0000003',    'BS_NAME'                          }; cnt = cnt+1;
      scan{cnt} = {'hC0000000',    'DATA_CHANNELS'                    }; cnt = cnt+1;
      scan{cnt} = {'hD0000000',    'DATA_CHANNEL'                     }; cnt = cnt+1;
      scan{cnt} = {'hD0000001',    'DATA_NAME'                        }; cnt = cnt+1;
      scan{cnt} = {'hD0000004',    'COLOR'	                          }; cnt = cnt+1;
      scan{cnt} = {'hD0000005',    'SAMPLETYPE'	                      }; cnt = cnt+1;
      scan{cnt} = {'hD0000006',    'BITS_PER_SAMPLE'                  }; cnt = cnt+1;
      scan{cnt} = {'hD0000007',    'RATIO_TYPE'                       }; cnt = cnt+1;
      scan{cnt} = {'hD0000008',    'RATIO_TRACK1'                     }; cnt = cnt+1;
      scan{cnt} = {'hD0000009',    'RATIO_TRACK2'                     }; cnt = cnt+1;
      scan{cnt} = {'hD000000A',    'RATIO_CHANNEL1'                   }; cnt = cnt+1;
      scan{cnt} = {'hD000000B',    'RATIO_CHANNEL2'                   }; cnt = cnt+1;
      scan{cnt} = {'hD000000C',    'RATIO_CONST1'                     }; cnt = cnt+1;
      scan{cnt} = {'hD000000D',    'RATIO_CONST2'                     }; cnt = cnt+1;
      scan{cnt} = {'hD000000E',    'RATIO_CONST3'                     }; cnt = cnt+1;
      scan{cnt} = {'hD000000F',    'RATIO_CONST4'                     }; cnt = cnt+1;
      scan{cnt} = {'hD0000010',    'RATIO_CONST5'                     }; cnt = cnt+1;
      scan{cnt} = {'hD0000011',    'RATIO_CONST6'                     }; cnt = cnt+1;
      scan{cnt} = {'hD0000012',    'RATIO_FIRST_IMAGES1'	            }; cnt = cnt+1;
      scan{cnt} = {'hD0000013',    'RATIO_FIRST_IMAGES2'	            }; cnt = cnt+1;
      scan{cnt} = {'hD0000014',    'DYE_NAME'                         }; cnt = cnt+1;
      scan{cnt} = {'hD0000015',    'DYE_FOLDER'                       }; cnt = cnt+1;
      scan{cnt} = {'hD0000016',    'SPECTRUM'                         }; cnt = cnt+1;
      scan{cnt} = {'hD0000017',    'ACQUIRE'                          }; cnt = cnt+1;
      scan{cnt} = {'h11000000',    'TIMERS'                           }; cnt = cnt+1;
      scan{cnt} = {'h12000000',    'TIMER'                            }; cnt = cnt+1;
      scan{cnt} = {'h12000001',    'TIMER_NAME'                       }; cnt = cnt+1;
      scan{cnt} = {'h12000003',    'TIMER_INTERVAL'                   }; cnt = cnt+1;
      scan{cnt} = {'h12000004',    'TIMER_TRIGGER_IN'                 }; cnt = cnt+1;
      scan{cnt} = {'h12000005',    'TIMER_TRIGGER_OUT'                }; cnt = cnt+1;
      scan{cnt} = {'h13000000',    'MARKERS'                          }; cnt = cnt+1;
      scan{cnt} = {'h14000000',    'MARKER'                           }; cnt = cnt+1;
      scan{cnt} = {'h14000001',    'MARKER_NAME'                      }; cnt = cnt+1;
      scan{cnt} = {'h14000002',    'MARKER_DESCRIPTION'               }; cnt = cnt+1;
      scan{cnt} = {'h14000003',    'MARKER_TRIGGER_IN'                }; cnt = cnt+1;
      scan{cnt} = {'h14000004',    'MARKER_TRIGGER_OUT'               }; cnt = cnt+1;
      scan{cnt} = {'hFFFFFFFF',    'END_SUBBLOCK'                     };
      obj.propertyList = scan;
    end
  end
  
end

