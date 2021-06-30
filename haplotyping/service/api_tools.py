from flask import Response, abort
from flask_restx import Namespace, Resource
import json, haplotyping

namespace = Namespace("tools", description="Several tools", path="/tools")

@namespace.route("/canonical/<kmer>")
@namespace.param("kmer", "k-mer")
class ToolsCanonicalSingle(Resource):
    @namespace.doc(description="Get canonical form for k-mer")
    def get(self, kmer):
        try:
            response = haplotyping.General.canonical(kmer)
            return Response(json.dumps(response), mimetype="application/json") 
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))

@namespace.route("/reverse-complement/<kmer>")
@namespace.param("kmer", "k-mer")
class ToolsReverseComplement(Resource):
    @namespace.doc(description="Get reverse-complement of k-mer")
    def get(self, kmer):
        try:
            response = haplotyping.General.reverse_complement(kmer)
            return Response(json.dumps(response), mimetype="application/json")
        except Exception as e:
            abort(e.code if hasattr(e,"code") else 500, str(e))
