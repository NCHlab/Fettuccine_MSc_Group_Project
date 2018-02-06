from flask import Blueprint, jsonify, request
from project import table_builder



tables = Blueprint('tables', __name__, url_prefix='/tables')

@tables.route("/erv_serverside_table", methods=['GET'])
def serverside_table_content_erv():
    data = table_builder.collect_data_serverside_erv(request)
    return jsonify(data)


@tables.route("/l1_serverside_table", methods=['GET'])
def serverside_table_content_l1():
    data = table_builder.collect_data_serverside_l1(request)
    return jsonify(data)
