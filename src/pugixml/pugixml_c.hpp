#include <cstring>

#include "pugixml.hpp"

using namespace pugi;

extern "C" {
  // xml_node functions

  const char* xml_node_name(xml_node_struct* node){
    return xml_node(node).name();
  }

  xml_node_struct* xml_node_child(xml_node_struct* node, char* name){
    return xml_node(node).child(name).internal_object();
  }

  xml_node_struct* xml_node_next_sibling(xml_node_struct* node, char* name){
    return xml_node(node).next_sibling(name).internal_object();
  }

  xml_attribute_struct* xml_node_attribute(xml_node_struct* node, char* name){
    return xml_node(node).attribute(name).internal_object();
  }

  const char* xml_node_child_value(xml_node_struct* node){
    return xml_node(node).child_value();
  }

  xml_node_struct* xml_node_text(xml_node_struct* node){
    return xml_node(node).internal_object();
  }

  // xml_attribute functions

  const char* xml_attribute_name(xml_attribute_struct* attribute){
    return xml_attribute(attribute).name();
  }

  const char* xml_attribute_value(xml_attribute_struct* attribute){
    return xml_attribute(attribute).value();
  }

  xml_attribute_struct* xml_attribute_next_attribute(xml_attribute_struct* attribute){
    return xml_attribute(attribute).next_attribute().internal_object();
  }

  bool xml_attribute_as_bool(xml_attribute_struct* attribute){
    return xml_attribute(attribute).as_bool();
  }

  int xml_attribute_as_int(xml_attribute_struct* attribute){
    return xml_attribute(attribute).as_int();
  }

  long long xml_attribute_as_llong(xml_attribute_struct* attribute){
    return xml_attribute(attribute).as_llong();
  }

  double xml_attribute_as_double(xml_attribute_struct* attribute){
    return xml_attribute(attribute).as_double();
  }

  // xml_text functions

  bool xml_text_as_bool(xml_node_struct* node){
    return xml_node(node).text().as_bool();
  }

  int xml_text_as_int(xml_node_struct* node){
    return xml_node(node).text().as_int();
  }

  long long xml_text_as_llong(xml_node_struct* node){
    return xml_node(node).text().as_llong();
  }

  double xml_text_as_double(xml_node_struct* node){
    return xml_node(node).text().as_double();
  }

  // xml_document functions

  xml_document* xml_document_load_file(char* filename){
    xml_document * doc = new xml_document();
    xml_parse_result result = doc->load_file(filename);
    return doc;
  }

  xml_node_struct* xml_document_document_element(xml_document* doc){
    return doc->document_element().internal_object();
  }

  void xml_document_clear(xml_document* doc){
    delete doc;
  }
}
