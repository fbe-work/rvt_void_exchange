"""
Commented 425 demo version of nestedDurchbruchCutter
Select a nested RVT link in project browser. The script will search through that
nested link. It will create accordingly sized void instances, based on the bounding
box size of existing "WD - RECH" generic models with parameter set to Leichtbauwand
from the nested RVT link. The Voids are set to the closest "Building Story" level.
It will then cut intersecting 10_AUS_WND_GKW/BST/DST walls with these voids.
prerequisite: ROHB_GEN_Aussparung_DD_rechteckig.rfa
source updated 20170508
"""

import clr
clr.AddReference("RevitAPI")
from Autodesk.Revit.DB import FilteredElementCollector as Fec
from Autodesk.Revit.DB import BuiltInCategory as Bic
from Autodesk.Revit.DB import Transaction, Outline, BoundingBoxIntersectsFilter, XYZ
from Autodesk.Revit.DB import InstanceVoidCutUtils, FamilySymbol, Structure, Options, BoundingBoxXYZ
from Autodesk.Revit.DB import ModelPathUtils, RevitLinkInstance, RevitLinkOptions
from Autodesk.Revit.DB import RevitLinkType
from Autodesk.Revit.DB.ModelPathUtils import ConvertModelPathToUserVisiblePath as ConvertToUserPath
from collections import defaultdict, namedtuple
import datetime
import sys
import os.path as op
from System.Diagnostics import Stopwatch
from rpw import doc, uidoc, db


def get_lvl_dicts(building_story=1):
    """
    returns dictionaries of levels, and levels with distance to above
    """
    LevelInfo = namedtuple("LevelInfo", "lvl_elev next_lvl_elev")
    lvl_dict = {}
    lvl_elevs = {}
    levels = Fec(doc).OfCategory(Bic.OST_Levels).WhereElementIsNotElementType().ToElements()
    for level in levels:
        if level.LookupParameter("Building Story").AsInteger() == building_story:
            # print("Level {0} found at elevation: {1}".format(level.Name, level.Elevation))
            lvl_elevs[level.Elevation] = level

    previous_lvl = None
    for elevation in sorted(lvl_elevs.keys(), reverse=True):
        current_lvl = lvl_elevs[elevation]
        if previous_lvl:
            elev_above_lvl = previous_lvl.Elevation
            lvl_dict[lvl_elevs[elevation].Name] = LevelInfo(elevation, elev_above_lvl)
            # print(lvl_elevs[elevation].Name, elevation, elev_above_lvl)
        else:
            lvl_dict[lvl_elevs[elevation].Name] = LevelInfo(elevation, 12)
            # print(lvl_elevs[elevation].Name, elevation, 12)
        previous_lvl = current_lvl

    return lvl_elevs, lvl_dict


def expand_contract_void(void, level_dict, expand=False):
    """
    expands contracts antennas of void according to level dict
    """
    void_top = void.Location.Point.Z
    lvl_below = get_closest_level_below(void.Location.Point, level_dict)
    lvl_below_elev = lvl_below.Elevation

    lvl_info = raw_next_lvl_elevations[level_dict[lvl_below_elev].Name]
    lvl_end = lvl_info.next_lvl_elev

    void_height = void.Symbol.LookupParameter("Aussp_Hoehe").AsDouble()
    threshold = 0.1
    top_line_length = lvl_end - void_top - threshold
    bottom_line_length = void_top - void_height - lvl_below_elev - threshold
    # print(lvl_start, lvl_end)
    # print(void_top, void_height)
    # print(top_line_length, bottom_line_length)
    if expand:
        if top_line_length > 0:
            void.LookupParameter("line_ext_top").Set(top_line_length)
        if bottom_line_length > 0:
            void.LookupParameter("line_ext_bottom").Set(bottom_line_length)
    else:
        void.LookupParameter("line_ext_top").Set(0.0)
        void.LookupParameter("line_ext_bottom").Set(0.0)


def timestamp_now():
    return datetime.datetime.now().strftime("%Y%m%d")


def get_bbox_of_solid_vertices(fam_inst):
    """
    retrieves a bounding box from first volumetric
    element in family instance
    """
    s0_geo = fam_inst.get_Geometry(geo_opt)

    pts = {"xmin": None, "xmax": None,
           "ymin": None, "ymax": None,
           "zmin": None, "zmax": None}

    for geo_elem in s0_geo:
        for solid in geo_elem.GetInstanceGeometry():
            if solid.GetType().Name == "Solid":
                if solid.Volume > 0.0:
                    edges = solid.Edges
                    for edge in edges:
                        ep = edge.AsCurve().GetEndPoint(0)
                        # first populate coords in dict
                        if not pts["xmin"]:
                            pts["xmin"] = ep.X

                        if not pts["xmax"]:
                            pts["xmax"] = ep.X

                        if not pts["ymin"]:
                            pts["ymin"] = ep.Y

                        if not pts["ymax"]:
                            pts["ymax"] = ep.Y

                        if not pts["zmin"]:
                            pts["zmin"] = ep.Z

                        if not pts["zmax"]:
                            pts["zmax"] = ep.Z

                        # then keep updating
                        if ep.X < pts["xmin"]:
                            pts["xmin"] = ep.X
                        elif ep.X > pts["xmax"]:
                            pts["xmax"] = ep.X

                        if ep.Y < pts["ymin"]:
                            pts["ymin"] = ep.Y
                        elif ep.Y > pts["ymax"]:
                            pts["ymax"] = ep.Y

                        if ep.Z < pts["zmin"]:
                            pts["zmin"] = ep.Z
                        elif ep.Z > pts["zmax"]:
                            pts["zmax"] = ep.Z
                    break

    vertex_bbox = BoundingBoxXYZ()
    vertex_bbox.Max = XYZ(pts["xmax"], pts["ymax"], pts["zmax"])
    vertex_bbox.Min = XYZ(pts["xmin"], pts["ymin"], pts["zmin"])
    return vertex_bbox


def mm_coord(coordinate):
    coord_x = int(coordinate.X * ft_mm)
    coord_y = int(coordinate.Y * ft_mm)
    coord_z = int(coordinate.Z * ft_mm)
    return XYZ(coord_x, coord_y, coord_z)


def get_closest_level_below(search_coord, lvl_elev_dict):
    """
    retrieves the closest level below a given coordinate
    """
    closest_elev_below = None
    for lvl_elev in sorted(lvl_elev_dict):
        if search_coord.Z > lvl_elev:
            closest_elev_below = lvl_elev
        else:
            break
    return lvl_elev_dict[closest_elev_below]


def get_local_matching_gen_mods(family_filter):
    """
    returns generic models filtered by linked source and family name
    """
    local_gen_mods = Fec(doc).OfCategory(Bic.OST_GenericModel).WhereElementIsNotElementType().ToElements()
    filtered_gen_mods = []
    for gen in local_gen_mods:
        if gen.LookupParameter("Linked_Source").AsString() == rvt_file_name:
            if gen.Symbol.FamilyName == family_filter:
                filtered_gen_mods.append(gen)
    return filtered_gen_mods


def correct_selection(selected_elems):
    """
    checks user selection
    """
    user_msg = "please select exactly one nested Revit link with Wanddurchbruch Instances"
    if selected_elems:
        if len(selected_elems) == 1:
            elem = selected_elems[0]
            if elem.IsNestedLink:
                if ".rvt" in elem.LookupParameter("Type Name").AsString():
                    print("single rvt nested link selected.")
                    return True
            else:
                print("selection was not a nested Link.")
                print(user_msg)
    else:
        print(user_msg)


def query_bbox_exists_in_cache(bbox_pt_tuple, bbox_cache):
    found = bbox_pt_tuple in bbox_cache
    # print("search in cache for: ", bbox_pt_tuple, found)
    return found


def bbx_ft_tuple_str(bbox):
    x_min = bbox.Min.X
    y_min = bbox.Min.Y
    z_min = bbox.Min.Z
    x_max = bbox.Max.X
    y_max = bbox.Max.Y
    z_max = bbox.Max.Z
    coord_tuple = (x_min, y_min, z_min, x_max, y_max, z_max)
    # print("ft cache: ", coord_tuple)
    return str(coord_tuple)


stopwatch = Stopwatch()
stopwatch.Start()

ignored_durchbruch_ids = [
    # "A001G01",
]

# general revitapi options and constants
ft_mm = 304.8
up = XYZ(0, 0, 1)
non_struct = Structure.StructuralType.NonStructural
geo_opt = Options()
rvt_link_opt = RevitLinkOptions(False)

selection = [doc.GetElement(elId) for elId in uidoc.Selection.GetElementIds()]

# create level overviews
lvl_elevations, next_lvl_elevations = get_lvl_dicts()
raw_lvl_elevations, raw_next_lvl_elevations = get_lvl_dicts(building_story=0)

# retrieve model elements
local_gen = Fec(doc).OfCategory(Bic.OST_GenericModel).WhereElementIsNotElementType().ToElements()
local_walls = Fec(doc).OfCategory(Bic.OST_Walls).WhereElementIsNotElementType().ToElements()
fam_symbols = Fec(doc).OfClass(FamilySymbol).ToElements()

# setup counter and variables
aussparung_types_dict = defaultdict(dict)
void_symbol = None
found_counter = 0
types_created = 0
existed_already_counter = 0
created_counter = 0
obsolete_counter = 0
family_filter_name = "ROHB_GEN_Aussparung_DD_rechteckig"
filtered_wall_names = ["10_AUS_WND_GK",
                       "10_AUS_BDB_BST",
                       "10_AUS_WND_DST",
                       ]

if correct_selection(selection):
    with db.Transaction("cur_nested_rvt_voids"):

        # create link instance of selected model link
        selected_link_type = selection[0]
        link_model_path = selected_link_type.GetExternalFileReference().GetAbsolutePath()
        print("found: {}".format(ConvertToUserPath(link_model_path)))
        rvt_link_type = RevitLinkType.Create(doc, link_model_path, rvt_link_opt)
        rvt_link = RevitLinkInstance.Create(doc, rvt_link_type.ElementId)
        rvt_doc = rvt_link.GetLinkDocument()
        rvt_file_name = op.split(rvt_doc.PathName)[-1]
        time_stamp = timestamp_now()

        # access local and linked generic models
        linked_gen_mods = Fec(rvt_doc).OfCategory(Bic.OST_GenericModel).WhereElementIsNotElementType().ToElements()
        existing_gen_mods = get_local_matching_gen_mods(family_filter_name)

        # setup caches
        bbx_ft_cache_existing = set()
        bbx_ft_cache_remote = set()

        # contrac void antennas
        for gen_mod in existing_gen_mods:
            expand_contract_void(gen_mod, raw_lvl_elevations)

        doc.Regenerate()

        print("creating cache")
        for gen_mod in existing_gen_mods:
            bbx = gen_mod.get_BoundingBox(None)
            bbx_ft_cache_existing.add(gen_mod.LookupParameter("Creating_bbox").AsString())
            print("added ", gen_mod.LookupParameter("Creating_bbox").AsString())

        # store the walls to be potentially cut
        filtered_wall_ids = []
        for wall in local_walls:
            if "WallType" in dir(wall):
                if wall.WallType.FamilyName != "Curtain Wall":
                    for filter_name in filtered_wall_names:
                        if wall.Name.startswith(filter_name):
                            filtered_wall_ids.append(wall.Id.IntegerValue)
        filtered_wall_ids = set(filtered_wall_ids)

        # store overview of existing void types
        for symbol in fam_symbols:
            if symbol.FamilyName == family_filter_name:
                symbol_name = symbol.LookupParameter("Type Name").AsString()
                void_symbol = symbol
                # print("Family Aussparung found:", symbol.FamilyName, symbol_name)
                aussparung_types_dict[symbol_name] = {"type_id": symbol.Id, "type_obj": symbol}

        # process all generic models from linked model
        for gen_mod in linked_gen_mods:

            # retrieve generic model parameter values
            fam_name = gen_mod.Symbol.FamilyName
            bauart_wand = gen_mod.LookupParameter("Bauart_Schott").AsString()
            durchbruch_nr = gen_mod.LookupParameter("Durchbruch Nummer").AsString()
            gewerk = gen_mod.LookupParameter("Gewerk").AsString()
            gen_lvl_name = rvt_doc.GetElement(gen_mod.LevelId).Name

            # filter generic models to be processed
            if not bauart_wand.endswith("ST"):
                continue

            if durchbruch_nr in ignored_durchbruch_ids:
                continue

            if not fam_name.startswith("WD - RECH"):
                continue

            print(70 * "_")
            # print(gen_mod)
            print(fam_name, gen_mod.Id.IntegerValue)

            # get bbox of edges from solid with volume
            void_bbox = get_bbox_of_solid_vertices(gen_mod)
            # print("found bbox: {} - {}".format(void_bbox.Min, void_bbox.Max))
            if not void_bbox:
                print("bbox not found!!")
                continue

            found_counter += 1
            coord = XYZ(void_bbox.Min.X, void_bbox.Min.Y, void_bbox.Max.Z)
            # print("mm_coord of void placement: {0}".format(mm_coord(coord)))

            # get closest level and offset to it
            closest_level = get_closest_level_below(coord, lvl_elevations)
            level_offset = coord - XYZ(0, 0, closest_level.Elevation)
            print("closest level: {}".format(closest_level.Name))
            # print("level offset: ", level_offset)

            # create coordinates check string
            # ask if bbox already exists
            remote_ft_bbx = bbx_ft_tuple_str(void_bbox)
            print("ask existing cache if bbx exists:", remote_ft_bbx)
            bbox_exists = query_bbox_exists_in_cache(remote_ft_bbx, bbx_ft_cache_existing)

            # add the bbox to remote bbox cache
            bbx_ft_cache_remote.add(remote_ft_bbx)

            # skip this void if it exists already in cache
            if bbox_exists:
                print("bbox exist already: {}".format(bbox_exists))
                existed_already_counter += 1
                continue

            # calculate void dimentsions
            ft_width = void_bbox.Max.X - void_bbox.Min.X
            ft_depth = void_bbox.Max.Y - void_bbox.Min.Y
            ft_height = void_bbox.Max.Z - void_bbox.Min.Z

            bb_width = ft_width * ft_mm
            bb_depth = ft_depth * ft_mm
            bb_height = ft_height * ft_mm

            name_width = str(int(round(bb_width)))
            name_depth = str(int(round(bb_depth)))
            name_height = str(int(round(bb_height)))

            # create according type name
            bb_target_type = name_width + "x" + name_depth + "x" + name_height
            print(bb_target_type,
                  "void width: " + name_width,
                  "void depth: " + name_depth,
                  "void height: " + name_height
                  )
            # print("bbox min coordinate: " + str(void_bbox.Min.ToString()))
            # print("bbox mm min coordinate: " + str(mm_coord(void_bbox.Min).ToString()))
            # print("bbox mm max coordinate: " + str(mm_coord(void_bbox.Max).ToString()))

            if "type_obj" in aussparung_types_dict[bb_target_type]:
                print("found type already")
                # print(aussparung_types_dict[bb_target_type]["type_obj"])
            elif not void_symbol:
                print("!!! please load the missing void cut family first - aborting. !!!")
                sys.exit()
            else:
                print("type not found! - creating: " + bb_target_type)
                dup_symbol = void_symbol.Duplicate(bb_target_type)
                aussparung_types_dict[bb_target_type] = {"type_id": dup_symbol.Id,
                                                         "type_obj": dup_symbol}
                dup_symbol.LookupParameter("Aussp_Breite").Set(ft_width)
                dup_symbol.LookupParameter("Aussp_Laenge").Set(ft_depth)
                dup_symbol.LookupParameter("Aussp_Hoehe").Set(ft_height)
                dup_symbol.Activate()
                doc.Regenerate()
                types_created += 1
                print("type created: " + bb_target_type)

            # create new void instance
            if not bbox_exists:
                if not aussparung_types_dict[bb_target_type]["type_obj"].IsActive:
                    aussparung_types_dict[bb_target_type]["type_obj"].Activate()

                gen_void = doc.Create.NewFamilyInstance(level_offset,
                                                        aussparung_types_dict[bb_target_type]["type_obj"],
                                                        closest_level,
                                                        non_struct)
                # mirror parameters
                gen_void.LookupParameter("Aussp_UK").Set(level_offset.Z - ft_height)
                gen_void.LookupParameter("TimeStamp").Set(time_stamp)
                gen_void.LookupParameter("Linked_Source").Set(rvt_file_name)
                gen_void.LookupParameter("Bauart_Schott").Set(bauart_wand)
                gen_void.LookupParameter("Durchbruch Nummer").Set(durchbruch_nr)
                gen_void.LookupParameter("Creating_bbox").Set(remote_ft_bbx)

                print("found params:", gewerk, bauart_wand, durchbruch_nr, gen_lvl_name,
                      "Id: ", gen_void.Id.IntegerValue)

                created_counter += 1

                # retrieve walls intersecting void bbox
                void_outline = Outline(void_bbox.Min, void_bbox.Max)
                bbox_filter = BoundingBoxIntersectsFilter(void_outline)
                cx_walls = Fec(doc).OfCategory(Bic.OST_Walls).WherePasses(bbox_filter).ToElements()

                doc.Regenerate()

                # cut intersecting walls with void
                for wall in cx_walls:
                    if wall.Id.IntegerValue in filtered_wall_ids:
                        print("attempt to cut the following elements: {} with {}".format(wall.Id,
                                                                                         gen_void.Id))
                        InstanceVoidCutUtils.AddInstanceVoidCut(doc, wall, gen_void)
                        gen_void.LookupParameter("Aussp_Typ").Set("WD")

        # mark obsolete voids with comment
        local_gen = get_local_matching_gen_mods(family_filter_name)
        for gen_mod in local_gen:
            bbx_ft_tuple = gen_mod.LookupParameter("Creating_bbox").AsString()
            bbox_not_obsolete = query_bbox_exists_in_cache(bbx_ft_tuple, bbx_ft_cache_remote)
            if not bbox_not_obsolete:
                gen_mod.LookupParameter("Comments").Set("void_obsolete: {}".format(time_stamp))
                print("not found -> obsolete: {}, coord: {}".format(gen_mod.Id.IntegerValue,
                                                                    bbx_ft_tuple))
                obsolete_counter += 1

            # expand void antennas
            expand_contract_void(gen_mod, raw_lvl_elevations, expand=True)

        # remove temporary link
        doc.Delete(rvt_link.Id)
        doc.Delete(rvt_link_type.ElementId)

# report for user
print(50 * "-")
print(str(types_created) + " new void size types created.")
print("{} voids found in linked model.".format(str(found_counter).zfill(4)))
print("{} voids were created.".format(str(created_counter).zfill(4)))
print("{} voids already existed.".format(str(existed_already_counter).zfill(4)))
print("{} voids are obsolete.".format(str(obsolete_counter).zfill(4)))
print("HdM_pyRevit nestedRevitDurchbruchCutter run in: ")

stopwatch.Stop()
timespan = stopwatch.Elapsed
print(timespan)

print("WARNING: Currently only creates orthogonal voids!!")
