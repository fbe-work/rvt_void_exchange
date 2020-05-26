"""
Select a nested RVT link. The script will search through this linked model.
It will create accordingly sized void instances, based on the bounding box size of
found "WD - RECH" generic models in that link, and will mirror a set of parameters.
The Voids are set to the closest "Building Story" level.
It will cut intersecting walls, floors, roofs and structural framings of predefined
types (see white lists and black lists) with these voids.
Furthermore will it check existing voids if they need to be recut, based on current
conditions in the model.
Rotation of void is currently not yet implemented!!
prerequisite: ROHB_GEN_Aussparung_DD_rechteckig.rfa
source: U:/Frederic_B/SCRIPTING/RVT/RPS/ROHB_GEN_Aussparung_DD_rechteckig.rfa
source updated 20200422
"""

import clr
clr.AddReference("RevitAPI")
from Autodesk.Revit.DB import FilteredElementCollector as Fec
from Autodesk.Revit.DB import FilteredWorksetCollector as Fwc
from Autodesk.Revit.DB import Outline, BoundingBoxIntersectsFilter, ElementIntersectsSolidFilter
from Autodesk.Revit.DB import XYZ, ElementId, Element
from Autodesk.Revit.DB import InstanceVoidCutUtils, Structure, Options, SolidOptions, BoundingBoxXYZ
from Autodesk.Revit.DB import ModelPathUtils, RevitLinkInstance, RevitLinkOptions, RevitLinkType, WorksetKind
from Autodesk.Revit.DB import Line, Curve, CurveLoop, GeometryCreationUtilities, ElementTransformUtils
from collections import defaultdict, namedtuple, deque
import datetime
import sys
import re
import time
import os
from System.Diagnostics import Stopwatch
from System.Collections.Generic import List
from rpw import doc, uidoc, db, DB
from pyrevit import script


# TODO check for lowest raw/non-raw level needed
# TODO implement rotation
# TODO check closest level via dict bounds


class LevelCatalogue:
    """
    Creates level catalogue with
    level sets specified by regex.
    Name and Elevation dicts are available
    for the sets and all collected levels.
    Levels are wrapped into namedtuple
    with additional attributes.

    if paired above or below
                          -> not raw_zone
    if raw zone above or below
      F-above / F-below   -> not valid
      F-above / R-below   -> raw_zone !!
      R-above / R-below   -> not valid
      R-above / F-below   -> not raw_zone
    """
    collected_level_elevations = None
    collected_levels_count = 0
    collected_levels = None
    paired_level_tolerance = 0.03
    raw_zones = []

    def __init__(self, name_regexes, levels=None):
        if not levels:
            levels = db.Collector(of_category="Levels", is_not_type=True).ToElements()
        desc_levels = sorted(levels, key=lambda x: x.Elevation, reverse=True)
        LevelInfo = namedtuple("LevelInfo", "Name Id lvl lvl_set Elevation Z lvl_above lvl_below paired_with")
        re_matched_levels = []
        all_collected_levels = []
        all_collected_lvl_ids = {}
        all_elevs_collector = {}
        lvl_pairs = {}

        lvl_sets = defaultdict(list)
        for lv in desc_levels:
            dprint(25 * "-")
            dprint(lv.Name)
            set_name, found_lvl = self._check_lvl_re(lv, name_regexes)
            if found_lvl:
                lvl_sets[set_name].append(found_lvl)
                re_matched_levels.append(found_lvl)

        for found_lvl in re_matched_levels:
            other_lvls = [lvl for lvl in re_matched_levels if lvl.Id.IntegerValue != found_lvl.Id.IntegerValue]
            pair_lvl = [level for level in other_lvls if
                        self._is_closer_than(self.paired_level_tolerance, level.Elevation, found_lvl.Elevation)]
            if pair_lvl:
                dprint("paired: ", found_lvl.Name, pair_lvl[0].Name)
                lvl_pairs[found_lvl.Id.IntegerValue] = pair_lvl[0].Id

        for set_name, levels in lvl_sets.items():
            levels.append(None)
            lvl_queue = deque([None], maxlen=3)
            level_collector = []
            name_collector = {}
            elevs_collector = {}

            print(35 * "-")
            print("LevelCatalogue, {}, {} levels".format(set_name, len(levels) - 1))
            for level in levels:
                # print("deque: ", len(lvl_queue), [l.Name for l in lvl_queue if "Name" in dir(l)])
                lvl_above, lvl_below = None, None
                lvl_queue.appendleft(level)

                paired_id = None

                if len(lvl_queue) < 1:
                    continue
                if not lvl_queue[1]:
                    continue

                lvl = lvl_queue[1]
                lvl_id = lvl.Id

                if lvl_pairs.get(lvl_id.IntegerValue):
                    paired_id = lvl_pairs.get(lvl_id.IntegerValue)

                lvl_below = lvl_queue[0]
                if len(lvl_queue) > 2:
                    lvl_above = lvl_queue[2]

                print(lvl.Name, lvl.Id.IntegerValue, lvl.Elevation)

                lvl_info = LevelInfo(lvl.Name, lvl_id, lvl, set_name,
                                     lvl.Elevation, lvl.Elevation, lvl_above,
                                     lvl_below, paired_id)

                level_collector.append(lvl_info)
                all_collected_levels.append(lvl_info)
                all_collected_lvl_ids[lvl.Id.IntegerValue] = lvl_info
                name_collector[lvl.Name] = lvl_info
                elevs_collector[lvl.Elevation] = lvl_info
                all_elevs_collector[lvl.Elevation] = lvl_info
                self.collected_levels_count += 1

            setattr(self, set_name, level_collector)
            setattr(self, set_name + "_names", name_collector)
            setattr(self, set_name + "_elevations", elevs_collector)

        self.collected_levels = sorted(all_collected_levels, key=lambda x: x.Elevation, reverse=True)
        self.collected_level_ids = all_elevs_collector
        self.collected_level_elevations = all_elevs_collector
        print(35 * "-")
        self._get_raw_zones()

    def _check_lvl_re(self, level, re_dicts):
        for re_dict in re_dicts:
            set_name = re_dict["name"]
            # dprint(set_name, re_dict["regex"].pattern)
            if re.match(re_dict["regex"], level.Name):
                dprint("{} matched: {}".format(re_dict["regex"].pattern, set_name))
                return set_name, level
        return None, None

    def _is_closer_than(self, tolerance, elev_a, elev_b):
        return abs(elev_a - elev_b) <= tolerance

    def _get_raw_zones(self):
        RawZone = namedtuple("RawZone", "name upper_bound lower_bound")
        for lvl in self.collected_levels:
            dprint(lvl.Elevation, lvl.Name)
            if getattr(lvl, "paired_with", None):
                dprint("paired_with: {}".format(lvl.paired_with))
                continue
            lvl_from_all_lvls_above = self.get_closest_lvl(lvl, self.collected_level_elevations, direction="above")
            if not lvl_from_all_lvls_above:
                continue
            above_set = lvl_from_all_lvls_above.lvl_set
            below_set = lvl.lvl_set
            above_elev = lvl_from_all_lvls_above.Elevation
            below_elev = lvl.Elevation
            dprint(lvl.Name, below_set, lvl_from_all_lvls_above.Name, above_set)
            if above_set == "OKFB" and below_set == "OKRB":
                zone_name = "{}::{}".format(lvl.Name, lvl_from_all_lvls_above.Name)
                print("raw_zone: {}".format(zone_name))
                self.raw_zones.append(RawZone(zone_name, above_elev, below_elev))
        print(35 * "-")

    def is_in_raw_zone(self, elevation):
        for zone in self.raw_zones:
            if zone.lower_bound < elevation < zone.upper_bound:
                return zone.name

    def get_closest_lvl(self, coordinate, lvl_elev_dict, direction=None):
        if direction == "above":
            closest_elev_below = None
            for lvl_elev in sorted(lvl_elev_dict, reverse=True):
                if coordinate.Z < lvl_elev:
                    closest_elev_below = lvl_elev
                else:
                    break
            return lvl_elev_dict.get(closest_elev_below)

        elif direction == "below":
            closest_elev_below = None
            # print("searching in: ", lvl_elev_dict)
            for lvl_elev in sorted(lvl_elev_dict):
                if coordinate.Z > lvl_elev:
                    closest_elev_below = lvl_elev
                else:
                    break
            return lvl_elev_dict.get(closest_elev_below)

        lvl_distances = {}
        for lvl_elev in lvl_elev_dict:
            distance = abs(coordinate.Z - lvl_elev)
            lvl_distances[distance] = lvl_elev_dict[lvl_elev]
        return lvl_distances[min(lvl_distances)]


def prerequisite_check():
    if not len(void_types):
        print("!!could not find void family: {}".format(family_filter_name))
        sys.exit()

    required_params = {
        'Linked_Source': {"category": 'Generic Models', "found": False},
        'Creating_bbox': {"category": 'Generic Models', "found": False},
    }
    model_params = {}
    param_binding_iter = doc.ParameterBindings.ForwardIterator()
    param_binding_iter.Reset()
    while param_binding_iter.MoveNext():
        param = param_binding_iter.Key
        binding = doc.ParameterBindings[param]
        categories = [cat.Name for cat in binding.Categories]
        model_params[param.Name] = categories
    for req_param_name, info in required_params.items():
        if model_params.get(req_param_name):
            # print("required_param {} exists".format(req_param_name))
            if info["category"] in model_params[req_param_name]:
                info["found"] = True
                # print("required param {} bound to {} found in project param bindings".format(req_param_name, info))
    for req_param_name, info in required_params.items():
        if not info["found"]:
            print("!!could not find param binding: {} for {}".format(req_param_name, info["category"]))


def expand_void_antenna(void, void_height, level_dict, expand_levels=None):
    threshold = 0.1
    if expand_levels:
        lvl_end        = expand_levels[0].Elevation
        lvl_below_elev = expand_levels[1].Elevation
        void_top = void.Location.Point.Z
        top_line_length = lvl_end - void_top - threshold
        bottom_line_length = void_top - void_height - lvl_below_elev - threshold
        # print("DD: {} - {}".format(top_line_length, bottom_line_length))
    else:
        void_top = void.Location.Point.Z
        lvl_below = lvl_cat.get_closest_lvl(void.Location.Point, level_dict, direction="below")
        lvl_below_elev = lvl_below.Elevation
        if 0.1 > void_top - lvl_below_elev > 0:
            # print("diff is: {}".format(void_top - lvl_below_elev))
            void_top = void.Location.Point.Z - 0.09
            lvl_below = lvl_cat.get_closest_lvl(void.Location.Point - XYZ(0, 0, -0.09), level_dict, direction="below")
            lvl_below_elev = lvl_below.Elevation

        lvl_end = lvl_below.lvl_above.Elevation
        top_line_length = lvl_end - void_top - threshold
        bottom_line_length = void_top - void_height - lvl_below_elev - threshold

    # dprint(lvl_below, lvl_end)
    # dprint(void_top, void_height)
    # dprint(top_line_length, bottom_line_length)
    if top_line_length > 0:
        void.LookupParameter("line_ext_top").Set(top_line_length)
    if bottom_line_length > 0:
        void.LookupParameter("line_ext_bottom").Set(bottom_line_length)


def timestamp_str_today():
    return datetime.datetime.now().strftime("%Y%m%d")


def get_bbox_of_solid_vertices(fam_inst):
    inst_geo = fam_inst.get_Geometry(geo_opt)
    ortho = True
    orthos = (
        XYZ.BasisX.ToString(), (XYZ.BasisX * -1).ToString(),
        XYZ.BasisY.ToString(), (XYZ.BasisY * -1).ToString(),
        XYZ.BasisZ.ToString(), (XYZ.BasisZ * -1).ToString(),
    )

    pts = {
        "xmin": None, "xmax": None,
        "ymin": None, "ymax": None,
        "zmin": None, "zmax": None,
    }

    for geo_elem in inst_geo:
        for solid in geo_elem.GetInstanceGeometry():
            if solid.GetType().Name == "Solid":
                if solid.Volume > 0.0:
                    faces = solid.Faces
                    for face in faces:
                        if ortho:
                            if face.FaceNormal.ToString() not in orthos:
                                ortho = False
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
                    if not ortho:
                        dprint(" non-ortho!! ")
                    break

    vertex_bbox = BoundingBoxXYZ()
    vertex_bbox.Max = XYZ(pts["xmax"], pts["ymax"], pts["zmax"])
    vertex_bbox.Min = XYZ(pts["xmin"], pts["ymin"], pts["zmin"])
    return vertex_bbox, ortho


def mm_coord(coordinate):
    coord_x = int(coordinate.X * ft_mm)
    coord_y = int(coordinate.Y * ft_mm)
    coord_z = int(coordinate.Z * ft_mm)
    return XYZ(coord_x, coord_y, coord_z)


def get_local_gen_mods_for_link(family_filter, link_file_name):
    local_gen_mods = db.Collector(of_category="GenericModel", is_not_type=True).ToElements()
    filtered_gen_mods = []
    for gen in local_gen_mods:
        if gen.LookupParameter("Linked_Source").AsString() == link_file_name:
            if gen.Symbol.FamilyName == family_filter:
                filtered_gen_mods.append(gen)
    return filtered_gen_mods


def get_filtered_element_instances(elements, filtered_element_names, blacklist=None):
    filtered_element_ids = List[ElementId]()
    for elem in elements:
        discard = False

        if elem.DesignOption:
            continue

        if "WallType" in dir(elem):
            if elem.WallType.FamilyName == "Curtain Wall":
                continue

        if blacklist:
            for fragment in blacklist:
                if fragment in elem.Name:
                    dprint(elem.Name, "elem discarded")
                    discard = True
            if discard:
                continue

        for filter_name in filtered_element_names:
            if elem.Name.startswith(filter_name):
                filtered_element_ids.Add(elem.Id)

    if filtered_element_ids:
        dprint("retrieved filtered elems: {} with {}".format(filtered_element_ids.Count, blacklist))
        return filtered_element_ids


def get_solid_extrusion_from_bbox(bbox):
    solid_opt = SolidOptions(ElementId.InvalidElementId, ElementId.InvalidElementId)
    b1 = XYZ(bbox.Min.X, bbox.Min.Y, bbox.Min.Z)
    b2 = XYZ(bbox.Max.X, bbox.Min.Y, bbox.Min.Z)
    b3 = XYZ(bbox.Max.X, bbox.Max.Y, bbox.Min.Z)
    b4 = XYZ(bbox.Min.X, bbox.Max.Y, bbox.Min.Z)
    bbox_height = bbox.Max.Z - bbox.Min.Z

    lines = List[Curve]()
    lines.Add(Line.CreateBound(b1, b2))
    lines.Add(Line.CreateBound(b2, b3))
    lines.Add(Line.CreateBound(b3, b4))
    lines.Add(Line.CreateBound(b4, b1))
    rectangle = [CurveLoop.Create(lines)]

    extrusion = GeometryCreationUtilities.CreateExtrusionGeometry(
        List[CurveLoop](rectangle),
        XYZ.BasisZ,
        bbox_height,
        solid_opt,
    )
    return extrusion


def verify_geometric_intersections(void_gen_mod_solid, cutting_already, potential_intersector_elems):
    if void_gen_mod_solid and potential_intersector_elems:
        cx_solid_filter = ElementIntersectsSolidFilter(void_gen_mod_solid)
        cx_ids = List[ElementId]()
        for elem in potential_intersector_elems:
            if not elem.Id.IntegerValue in cutting_already:
                cx_ids.Add(elem.Id)
        if cx_ids:
            geo_cx = Fec(doc, cx_ids).WherePasses(cx_solid_filter).ToElements()
            return geo_cx
    return List[Element]()


def cut_solid_cx_elems(solid_cx_elems, cat, void, schlitz, cut_types):
    dprint((solid_cx_elems, cat, void, schlitz, cut_types))
    void_cuts_types = None
    void_cuts = ""
    for elem in solid_cx_elems:
        dprint(elem.Id.IntegerValue)
        print("attempt to cut the following elements: "
              "{}_id: {} with void_id: {}".format(
            cat,
            output.linkify(elem.Id),
            output.linkify(void.Id))
        )
        InstanceVoidCutUtils.AddInstanceVoidCut(doc, elem, void)
        void_cuts = cat
        void_cuts_types = cut_types[schlitz][cat]

    return void_cuts, void_cuts_types


def get_solid_intersecting_elems(cut_categories, filtered_elem_ids, bbox_filter):
    solid_cut_elems = {}
    for cat in cut_categories:
        if filtered_elem_ids[cat]:
            cx_elems = Fec(doc, filtered_elem_ids[cat]).WherePasses(bbox_filter).ToElements()
            solid_cut_elems[cat] = verify_geometric_intersections(void_solid, [], cx_elems) or []
            dprint(cat, [elem.Id.IntegerValue for elem in cx_elems])
        else:
            solid_cut_elems[cat] = []
    return solid_cut_elems


def uncut_elems_cut_by_void(void):
    dprint("void {} uncutting.".format(void.Id.IntegerValue))
    cut_elem_ids = InstanceVoidCutUtils.GetElementsBeingCut(void)
    if cut_elem_ids:
        for cut_elem_id in cut_elem_ids:
            cut_elem = doc.GetElement(cut_elem_id)
            InstanceVoidCutUtils.RemoveInstanceVoidCut(doc, cut_elem, void)


def set_void_viz_params(fam_inst, mode=None):
    dprint(mode)
    set_param = {
        "DD": "VIS_Bodendurchbruch",
        "BS": "VIS_Bodenschlitz",
        "WD": "VIS_Wanddurchbruch",
        "WS": "VIS_Wandschlitz",
        "TD": "VIS_Wanddurchbruch",
        "TS": "VIS_Wandschlitz",
    }

    param_set = set(set_param.values())
    param_set.remove(set_param[mode])
    dprint("excluded from switching off:", mode)

    if set_param.get(mode):
        for param_name in param_set:
            fam_inst.LookupParameter(param_name).Set(0)


def correct_selection(selected_elems):
    user_msg = "aborted: please select exactly one nested Revit link with Wanddurchbruch Instances."
    if selected_elems:
        if len(selected_elems) == 1:
            elem = selected_elems[0]
            if elem.Category.Name != "RVT Links":
                print(user_msg)
                return
            elif elem.IsNestedLink:
                if ".rvt" in elem.LookupParameter("Type Name").AsString():
                    print("single rvt nested link selected.")
                    return True
            else:
                print("selection was not a nested Link.")
                print(user_msg)
    else:
        print(user_msg)


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


def get_tuple_from_bbx_str(coord_str):
    return tuple([float(num) for num in coord_str[1:-1].split(",")])


def bbox_centroid(bbox):
    """
    Retrieves centroid of a bbox.
    :param bbox:
    :return: XYZ(): centroid point
    """
    centroid = XYZ((bbox.Max.X - bbox.Min.X) / 2 + bbox.Min.X,
                         (bbox.Max.Y - bbox.Min.Y) / 2 + bbox.Min.Y,
                         (bbox.Max.Z - bbox.Min.Z) / 2 + bbox.Min.Z)
    return centroid


def dprint(*args, **kwargs):
    """
    prints args and kwargs only if the script
    cas called with debug ( ctrl + click )
    :param args:
    :param kwargs:
    :return:
    """
    # noinspection PyUnresolvedReferences
    if __forceddebugmode__:
        print(args, kwargs)


def get_param_val(elem, param_name, param=None):
    """
    Retrieves parameter value of element or parameter
    or its standard empty value for its type.
    :param elem: the element holding the parameter
    :param param_name: name of the parameter
    :param param: optionally the param instead of elem
    :return: value of the parameter or empty of type
    """
    if not param:
        param = elem.LookupParameter(param_name)
    if param:
        dtype = param.StorageType
        if param.HasValue:
            return dtype_methods[dtype](param)
        return dtype_empty[dtype]


stopwatch = Stopwatch()
stopwatch.Start()

dtype_methods = {
    DB.StorageType.String   : DB.Parameter.AsString,
    DB.StorageType.Integer  : DB.Parameter.AsDouble,
    DB.StorageType.Double   : DB.Parameter.AsDouble,
    DB.StorageType.ElementId: DB.Parameter.AsElementId,
}
dtype_empty = {
    DB.StorageType.String : "",
    DB.StorageType.Integer: 0.0,
    DB.StorageType.Double : 0.0,
    DB.StorageType.ElementId: DB.ElementId(-1),
}

ignored_durchbruch_ids = [
    # "S037U01",  # UG1
]
ft_mm = 304.8
up = XYZ(0, 0, 1)
non_struct = Structure.StructuralType.NonStructural
geo_opt = Options()
rvt_link_opt = RevitLinkOptions(False)

project_lvl_code = None
project_number = doc.ProjectInformation.LookupParameter("Project Number").AsString()
if project_number.isdigit() and len(project_number) == 3:
    project_lvl_code = "B{project_number}".format(project_number=project_number[-2:])
print("project: {project_lvl_code}, level prefix is: {project_lvl_code}".format(project_number=project_number,
                                                                                project_lvl_code=project_lvl_code))

re_levels = [
    {
        "name": "OKRB",
        "regex": re.compile("{}.+OKRB.?".format(project_lvl_code)),
    },
    {
        "name": "OKFB",
        "regex": re.compile("{}".format(project_lvl_code)),
    },
]

selection = [doc.GetElement(elId) for elId in uidoc.Selection.GetElementIds()]
lvl_cat = LevelCatalogue(re_levels)

ws_name_map = {ws.Name: ws for ws in Fwc(doc).OfKind(WorksetKind.UserWorkset).ToWorksets()}
void_ws = ws_name_map.get("F90_SD")

family_filter_name = "ROHB_GEN_Aussparung_DD_rechteckig"

local_walls        = db.Collector(of_category="Walls",             is_not_type=True ).ToElements()
local_floors       = db.Collector(of_category="Floors",            is_not_type=True ).ToElements()
local_roofs        = db.Collector(of_category="Roofs",             is_not_type=True ).ToElements()
local_str_framings = db.Collector(of_category="StructuralFraming", is_not_type=True ).ToElements()
void_types         = db.Collector(of_category="GenericModel",      is_type=True,
                                  where=lambda x: x.FamilyName == family_filter_name).ToElements()

counter = defaultdict(int)
aussparung_types_dict = defaultdict(dict)
local_void_by_bbx = {}
void_not_cut_ids = []
void_symbol = None

cut_categories = ["walls", "floors", "roofs", "structural_framing"]

allowed_types = [
    "GK",
    "BS",
    "DS",
    "FC",
]
filtered_wall_names = [
    "XX_",
    "10_",
    "10_",
    "20_",
]
filtered_roof_names = [
    "10_",
]
filtered_str_framings_names = [
    "XX_",
]
filtered_floor_names = [
    "XX_",
    "10_",
]
floor_names_blacklist = [
    "_HOB_",
    "_DOB_",
]

discipline_void_fam_re_pat = {
    "F50": re.compile(r"50_AUS_ELT_[WD]D - v.+"),
    "F90": re.compile(r"(WD|DD|WS|BS) - RECH \(whp v6a\)"),
}

void_cut_types = {
    False: {
        "walls": "WD",
        "floors": "DD",
        "roofs": "DD",
        "structural_framing": "TD",
    },
    True: {
        "walls": "WS",
        "floors": "BS",
        "roofs": "BS",
        "structural_framing": "TS",
    },
}
void_cut_type_names = [
    "WD",
    "DD",
    "WS",
    "BS",
    "TD",
    "TS",
]

collect_void_params = {
    "Durchbruch Nummer",
    "Gewerk",
    "Verortung",
    "SD_Verortung",
}

for void_type in void_types:
    symbol_name = void_type.LookupParameter("Type Name").AsString()
    # dprint("Family Aussparung found:", void_type.FamilyName, symbol_name)
    aussparung_types_dict[symbol_name] = {"type_id": void_type.Id, "type_obj": void_type}
    void_symbol = void_type

void_bbx_map    = {}
void_height_map = {}

output = script.get_output()

dprint("model check:")
prerequisite_check()

if not correct_selection(selection):
    sys.exit()

with db.Transaction("cut_nested_rvt_voids"):

    selected_link_type = selection[0]
    link_model_path = selected_link_type.GetExternalFileReference().GetAbsolutePath()
    link_user_path = ModelPathUtils.ConvertModelPathToUserVisiblePath(link_model_path)

    discipline_void_re_pat = ""

    for discipline in discipline_void_fam_re_pat:
        if discipline in link_user_path:
            discipline_void_re_pat = discipline_void_fam_re_pat[discipline]
            break

    print("found: {}".format(link_user_path))
    print("date modified: {}".format(time.ctime(os.path.getmtime(link_user_path))))

    rvt_link_type = RevitLinkType.Create(doc, link_model_path, rvt_link_opt)
    link_inst = RevitLinkInstance.Create(doc, rvt_link_type.ElementId)

    rvt_link = link_inst
    rvt_doc = rvt_link.GetLinkDocument()
    linked_gen_mods = db.Collector(doc=rvt_doc, of_category="GenericModel", is_not_type=True).ToElements()
    rvt_file_name = os.path.split(rvt_doc.PathName)[-1]
    time_stamp = timestamp_str_today()

    existing_gen_mods = get_local_gen_mods_for_link(family_filter_name, rvt_file_name)

    bbx_ft_cache_existing = set()
    bbx_ft_cache_remote   = set()

    print("creating local voids cache...")
    for gen_mod in existing_gen_mods:
        existing_bbx_str = gen_mod.LookupParameter("Creating_bbox").AsString()
        bbx_tuple = get_tuple_from_bbx_str(existing_bbx_str)
        bbx_ft_cache_existing.add(existing_bbx_str)
        dprint("added ", existing_bbx_str)
        local_void_by_bbx[existing_bbx_str] = gen_mod
        void_height_map[gen_mod.Id.IntegerValue] = bbx_tuple[-1] - bbx_tuple[-4]

    print("creating remote voids cache...")
    for gen_mod in linked_gen_mods:
        fam_name = gen_mod.Symbol.FamilyName
        if not re.match(discipline_void_re_pat, fam_name):
            continue
        void_bbox, void_ortho = get_bbox_of_solid_vertices(gen_mod)
        if void_bbox:
            dprint("remote cache id:", gen_mod.Id.IntegerValue)
            void_bbx_map[gen_mod.Id.IntegerValue] = {"bbox": void_bbox, "ortho": void_ortho}
            remote_ft_bbx = bbx_ft_tuple_str(void_bbox)
            bbx_ft_cache_remote.add(remote_ft_bbx)

    obsolete_coords = bbx_ft_cache_existing - bbx_ft_cache_remote
    obsolete_ids = set()
    for coord in obsolete_coords:
        counter["voids_existing_obsolete"] += 1
        dprint(coord)
        obsolete_void = local_void_by_bbx[coord]
        obsolete_ids.add(obsolete_void.Id)
        obsolete_void.LookupParameter("Comments").Set("void_obsolete: {}".format(time_stamp))
        uncut_elems_cut_by_void(obsolete_void)

    filtered_elem_ids = {
        "walls":              get_filtered_element_instances(local_walls, filtered_wall_names),
        "floors":             get_filtered_element_instances(local_floors, filtered_floor_names,
                                                             blacklist=floor_names_blacklist),
        "roofs":              get_filtered_element_instances(local_roofs, filtered_roof_names),
        "structural_framing": get_filtered_element_instances(local_str_framings, filtered_str_framings_names),
    }

    prog_bar_last_step = int(0.1 * linked_gen_mods.Count)
    prog_bar_total = linked_gen_mods.Count + prog_bar_last_step
    i = 0  # instead exit out if linked_gen_mods.Count == 0

    for i, gen_mod in enumerate(linked_gen_mods, 1):
        output.update_progress(i, prog_bar_total)
        fam_name = gen_mod.Symbol.FamilyName

        host_info = ""
        is_schlitz = False

        if not re.match(discipline_void_re_pat, fam_name):
            continue

        print(70 * "_")
        dprint(fam_name, gen_mod.Id.IntegerValue)

        durchbruch_nr = gen_mod.LookupParameter("Durchbruch Nummer").AsString()
        if durchbruch_nr in ignored_durchbruch_ids:
            continue
        print("rvt id: {} , Durchbruch Nummer: {}".format(gen_mod.Id.IntegerValue, durchbruch_nr))

        gewerk = gen_mod.LookupParameter("Gewerk").AsString()

        verortung = gen_mod.LookupParameter("Verortung")
        if verortung:
            verortung = verortung.AsString()

        sd_verortung = gen_mod.LookupParameter("SD_Verortung")
        if sd_verortung:
            sd_verortung = sd_verortung.AsString()

        host_info = verortung or sd_verortung or ""
        dprint(host_info)

        if discipline == "F50" and sd_verortung:
            if "schlitz" in sd_verortung:
                is_schlitz = True
        elif discipline == "F90":
            if fam_name.startswith("WS") or fam_name.startswith("BS"):
                is_schlitz = True

        dprint("is_schlitz: {}".format(is_schlitz), discipline)
        if is_schlitz:
            print("is_schlitz: {}".format(is_schlitz), discipline)

        void_bbox  = void_bbx_map[gen_mod.Id.IntegerValue]["bbox"]
        void_ortho = void_bbx_map[gen_mod.Id.IntegerValue]["ortho"]
        dprint("found bbox: {} - {}".format(void_bbox.Min, void_bbox.Max))

        if not void_bbox:
            print("bbox not found!!")

        if void_bbox:
            gen_void = None
            recut_existing = None
            cx_cut_elems = {category: [] for category in cut_categories}
            counter["voids_linked_found"] += 1
            coord = XYZ(void_bbox.Min.X, void_bbox.Min.Y, void_bbox.Max.Z)
            dprint("mm_coord of void placement: {0}".format(mm_coord(coord)))

            closest_lvl_below = lvl_cat.get_closest_lvl(coord, lvl_cat.OKFB_elevations, direction="below")
            placement_level = closest_lvl_below
            level_offset = coord - XYZ(0, 0, closest_lvl_below.Elevation)
            dprint("closest level: {}".format(closest_lvl_below.Name))
            dprint("level offset: ", level_offset)
            is_raw_zone = lvl_cat.is_in_raw_zone(void_bbox.Max.Z)

            remote_ft_bbx = bbx_ft_tuple_str(void_bbox)
            dprint("ask existing cache if bbx exists:", remote_ft_bbx)
            bbox_exists = remote_ft_bbx in bbx_ft_cache_existing
            bbx_ft_cache_remote.add(remote_ft_bbx)

            void_cuts = {category: False for category in cut_categories}
            void_cut_type = ""

            void_outline = Outline(void_bbox.Min, void_bbox.Max)
            void_solid   = get_solid_extrusion_from_bbox(void_bbox)
            bbox_filter  = BoundingBoxIntersectsFilter(void_outline, -0.03)

            if bbox_exists:
                dprint("bbox exist already: {}".format(bbox_exists))
                counter["voids_existing_already"] += 1
                local_existing_void = local_void_by_bbx[remote_ft_bbx]
                local_existing_void_id = local_existing_void.Id.IntegerValue
                local_existing_void.LookupParameter("Durchbruch Nummer").Set(durchbruch_nr)
                print("bbox exist already, updating void data, checking for recut.")

                already_cutting_ids = {eid.IntegerValue for eid in InstanceVoidCutUtils.GetElementsBeingCut(local_existing_void)}

                if not already_cutting_ids:
                    print(" -> existing void {} cuts nothing currently!!".format(local_existing_void_id))
                    counter["voids_existing_cuts_nothing"] += 1

                for cat in cx_cut_elems:
                    if filtered_elem_ids[cat]:
                        cx_elems = Fec(doc, filtered_elem_ids[cat]).WherePasses(bbox_filter).ToElements()
                        cx_cut_elems[cat] = verify_geometric_intersections(void_solid, already_cutting_ids,
                                                                           cx_elems) or []

                for cat in cx_cut_elems:
                    for elem in cx_cut_elems[cat]:
                        cut_elem_id = elem.Id.IntegerValue
                        note = ""
                        if local_existing_void_id < cut_elem_id:
                            note = ", elem to cut younger than void!!"
                            counter["voids_existing_recut_elem_younger"] += 1
                        else:
                            counter["voids_existing_recut"] += 1
                        print(" void: {} not cutting: {}{}".format(local_existing_void_id, cut_elem_id, note))
                        recut_existing = local_existing_void
                        dprint("attempt to recut: {}".format(recut_existing))
                        # mixed cut/uncut state so we reset to all not cur and recut
                        uncut_elems_cut_by_void(recut_existing)
                        # cut state has changed now, so we need to re-evaluate solid intersects right before recut

            ft_width = void_bbox.Max.X - void_bbox.Min.X
            ft_depth = void_bbox.Max.Y - void_bbox.Min.Y
            ft_height = void_bbox.Max.Z - void_bbox.Min.Z

            bb_width = ft_width * ft_mm
            bb_depth = ft_depth * ft_mm
            bb_height = ft_height * ft_mm

            name_width  = str(int(round(bb_width)))
            name_depth  = str(int(round(bb_depth)))
            name_height = str(int(round(bb_height)))

            bb_target_type = name_width + "x" + name_depth + "x" + name_height

            dprint("bbox min coordinate: "    + str(         void_bbox.Min.ToString()))
            dprint("bbox mm min coordinate: " + str(mm_coord(void_bbox.Min).ToString()))
            dprint("bbox mm max coordinate: " + str(mm_coord(void_bbox.Max).ToString()))

            if "type_obj" in aussparung_types_dict[bb_target_type]:
                dprint("found type already")
                dprint(aussparung_types_dict[bb_target_type]["type_obj"])
            else:
                print("type not yet existing, creating: " + bb_target_type)
                dup_symbol = void_symbol.Duplicate(bb_target_type)
                aussparung_types_dict[bb_target_type] = {"type_id": dup_symbol.Id,
                                                         "type_obj": dup_symbol}
                dup_symbol.LookupParameter("Aussp_Breite").Set(ft_width)
                dup_symbol.LookupParameter("Aussp_Laenge").Set(ft_depth)
                dup_symbol.LookupParameter("Aussp_Hoehe" ).Set(ft_height)
                dup_symbol.Activate()
                doc.Regenerate()
                counter["void_types_created"] += 1
                print("type created: " + bb_target_type)

            if not bbox_exists:
                if not aussparung_types_dict[bb_target_type]["type_obj"].IsActive:
                    aussparung_types_dict[bb_target_type]["type_obj"].Activate()

                gen_void = doc.Create.NewFamilyInstance(
                    level_offset,
                    aussparung_types_dict[bb_target_type]["type_obj"],
                    placement_level.lvl,
                    non_struct,
                )

                void_height_map[gen_void.Id.IntegerValue] = ft_height

                gen_void.LookupParameter("TimeStamp"        ).Set(time_stamp)
                gen_void.LookupParameter("Linked_Source"    ).Set(rvt_file_name)
                gen_void.LookupParameter("Creating_bbox"    ).Set(remote_ft_bbx)
                gen_void.LookupParameter("Durchbruch Nummer").Set(durchbruch_nr)
                gen_void.LookupParameter("Aussp_UK"         ).Set(level_offset.Z - ft_height)
                gen_void.LookupParameter("Aussp_Host"       ).Set(host_info)
                # gen_void.LookupParameter("Bauart_Schott"       ).Set(bauart_wand)
                if is_raw_zone:
                    level_above = lvl_cat.get_closest_lvl(coord, lvl_cat.OKFB_elevations, direction="above")
                    above_level_offset = coord - XYZ(0, 0, level_above.Elevation)
                    gen_void.LookupParameter("Aussp_UK").Set(above_level_offset.Z - ft_height)
                    print("raw offset: {}".format((above_level_offset.Z - ft_height) * ft_mm))

                counter["voids_new_created"] += 1
                if not void_ortho:
                    found_rotation = gen_mod.Location.Rotation
                    print("void is not ortho, angle: {}".format(found_rotation))

                if not void_ortho:
                    loc_pt = bbox_centroid(void_bbox)
                    rotation_angle = gen_mod.Location.Rotation
                    # zax = Line.CreateBound(loc_pt, loc_pt + XYZ(0, 0, 1))
                    # ElementTransformUtils.RotateElement(doc, gen_void.Id, zax, rotation_angle)

                if void_ws:
                    gen_void.LookupParameter("Workset").Set(void_ws.Id.IntegerValue)

                doc.Regenerate()

            if recut_existing:
                gen_void = recut_existing
                print("recutting: {}".format(recut_existing.Id.IntegerValue))

            if gen_void:
                cx_cut_elems = get_solid_intersecting_elems(cut_categories, filtered_elem_ids, bbox_filter)
                for cat in cx_cut_elems:
                    void_cut_cat, void_cut_type_found = cut_solid_cx_elems(
                        cx_cut_elems[cat],
                        cat,
                        gen_void,
                        is_schlitz,
                        void_cut_types,
                    )
                    void_cuts[void_cut_cat] = True
                    if void_cut_type_found:
                        void_cut_type += void_cut_type_found + " "
                        dprint(" !!found: ", void_cuts, void_cut_type)

                # void cut info
                if len(void_cut_type.strip()) == 2:
                    void_cut_type = void_cut_type.strip()

                if void_cuts["walls"]:
                    dprint("apply Aussp_Host_Type info")
                    re_wall_ausbau = re.compile("10_AUS_WND_(.+)")
                    wall_cut_types= []
                    found_at_wall_type = None
                    for cut_wall in cx_cut_elems["walls"]:
                        wall_type_name = cut_wall.Name
                        dprint("Aussp_Host_Type info: {}".format(wall_type_name))
                        re_found = re.findall(re_wall_ausbau, wall_type_name)
                        if re_found:
                            found_at_wall_type = re_found
                    if found_at_wall_type:
                        dprint("writing Aussp_Host_Type info: {}".format(found_at_wall_type[0]))
                        gen_void.LookupParameter("Aussp_Host_Type").Set(found_at_wall_type[0])

                if void_cuts["walls"] and void_cuts["floors"]:
                    void_cut_type = "WD/WS/DD/BS"
                    print("!! THIS VOID CUTS WALLS AND FLOORS/ROOFS !!")
                    gen_void.LookupParameter("Aussp_Typ"  ).Set(void_cut_type)
                    gen_void.LookupParameter("Aussp_Debug").Set(
                        "{}: did cut multiple categories!! engineer_info on this:{}".format(time_stamp, host_info)
                    )

                elif set(void_cuts.values()) == {False}:
                    void_cut_type = "no_cut!"
                    print("!! THIS VOID DOES NOT CUT ANYTHING !!")
                    gen_void.LookupParameter("Aussp_Debug").Set(
                            "{}: did not cut anything!! engineer_info on this:{}".format(time_stamp, host_info)
                    )
                    counter["voids_new_cuts_nothing"] += 1
                    void_not_cut_ids.append(gen_void.Id)

                if void_cut_type in void_cut_type_names:
                    dprint(void_cut_type)
                    gen_void.LookupParameter("Aussp_Typ").Set(void_cut_type)
                    set_void_viz_params(gen_void, void_cut_type)

                debug_param_val = get_param_val(gen_void, "Aussp_Debug")
                if not void_ortho and "non-ortho" not in debug_param_val:
                    gen_void.LookupParameter("Aussp_Debug").Set("non-ortho " + debug_param_val)

    # antenna expand
    local_gen = get_local_gen_mods_for_link(family_filter_name, rvt_file_name)
    for gen_mod in local_gen:
        void_cuts_floors = get_param_val(gen_mod, "Aussp_Typ") == "DD" or \
                           get_param_val(gen_mod, "Aussp_Typ") == "BS"
        void_height = void_height_map[gen_mod.Id.IntegerValue]
        if void_cuts_floors:
            closest_raw_lvl = lvl_cat.get_closest_lvl(gen_mod.Location.Point, lvl_cat.OKRB_elevations)
            closest_raw_lvl_z_coord = XYZ(0, 0, closest_raw_lvl.Elevation)
            raw_lvl_above = lvl_cat.get_closest_lvl(closest_raw_lvl_z_coord, lvl_cat.OKRB_elevations, direction="above")
            raw_lvl_below = lvl_cat.get_closest_lvl(closest_raw_lvl_z_coord, lvl_cat.OKRB_elevations, direction="below")

            # print("expand to: ", closest_raw_level.Elevation)
            expand_void_antenna(gen_mod, void_height, lvl_cat.OKRB_elevations,
                                expand_levels=(raw_lvl_above, raw_lvl_below))

            dprint("-found closest raw levels: {}\n{}\n{}".format(
                closest_raw_lvl.Name,
                raw_lvl_above.Name,
                raw_lvl_below.Name,
            ))

            gen_mod.LookupParameter("Closest_raw_level").Set(closest_raw_lvl.Name)
        else:
            expand_void_antenna(gen_mod, void_height, lvl_cat.OKRB_elevations)

    doc.Delete(link_inst.Id)
    doc.Delete(rvt_link_type.ElementId)
    output.update_progress(prog_bar_total, prog_bar_total)

dprint(70 * "-")
dprint("{} coordinates of obsolete voids:".format(counter["voids_existing_obsolete"]))
for coord in obsolete_coords:
    dprint(coord)

print(70 * "-")
print("ids of {} obsolete voids:".format(counter["voids_existing_obsolete"]))
for void_id in obsolete_ids:
    print(output.linkify(void_id))
print(70 * "-")
print("new voids that did not cut anything:")
print([output.linkify(elem_id) for elem_id in void_not_cut_ids])

print(70 * "-")
for topic in sorted(counter):
    print("{} : {}".format(str(counter[topic]).zfill(4), topic))

print("\nHdM_pyRevit nestedRevitDurchbruchCutter run in: ")
stopwatch.Stop()
timespan = stopwatch.Elapsed
print(timespan)

print("WARNING: This script currently only creates orthogonal voids!!")
