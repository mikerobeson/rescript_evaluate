# ----------------------------------------------------------------------------
# Copyright (c) 2019-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import importlib

from q2_types.genome_data import GenomeData, Loci, Proteins, Genes, DNASequence
from q2_types.metadata import ImmutableMetadata
from qiime2.core.type import TypeMatch
from qiime2.plugin import (Str, Plugin, Choices, List, Citations, Range, Int,
                           Float, Visualization, Bool, TypeMap, Metadata,
                           MetadataColumn, Categorical, Numeric, Collection)

from .cross_validate import (evaluate_cross_validate,
                             evaluate_fit_classifier)

from q2_types.feature_data import (FeatureData, Taxonomy, Sequence,
                                    AlignedSequence, RNASequence,
                                    AlignedRNASequence, ProteinSequence)
from q2_types.tree import Phylogeny, Rooted
from q2_feature_classifier.classifier import (_parameter_descriptions,
                                              _classify_parameters)
from q2_feature_classifier._taxonomic_classifier import TaxonomicClassifier
import rescript_eval
import rescript

from rescript.types._type import SILVATaxonomy, SILVATaxidMap
from rescript.types.methods import reverse_transcribe
from rescript.ncbi import (
    get_ncbi_data, _default_ranks, _allowed_ranks, get_ncbi_data_protein)


citations = Citations.load('citations.bib', package='rescript_eval')

plugin = Plugin(
    name='rescript_eval',
    version=rescript.__version__,
    website="https://github.com/nbokulich/RESCRIPt",
    package='rescript',
    description=('Supplementary RESCRIPt tools for '
                 'reference sequence and taxonomy evaluation.'),
    short_description=(
        'Supplementary RESRIPt evaluation tools.'),
    citations=[citations['Robeson2021rescript']]
)


LINEPLOT_XAXIS_INTERPRETATION = (
    'The x-axis in these plots represents the taxonomic '
    'levels present in the input taxonomies so are labeled numerically '
    'instead of by rank, but typically for 7-level taxonomies these will '
    'represent: 1 = domain/kingdom, 2 = phylum, 3 = class, 4 = order, '
    '5 = family, 6 = genus, 7 = species.')

rank_handle_description = (
    'Regular expression indicating which taxonomic rank a label '
    'belongs to; this handle is stripped from the label '
    'prior to operating on the taxonomy. The net '
    'effect is that ambiguous or empty levels can be '
    'removed prior to comparison, enabling selection of '
    'taxonomies with more complete taxonomic information. '
    'For example, "^[dkpcofgs]__" will recognize greengenes or silva rank '
    'handles. ')

rank_handle_extra_note = (
    'Note that rank_handles are removed but not replaced; use the '
    'new_rank_handle parameter to replace the rank handles.')

labels_description = (
    'List of labels to use for labeling evaluation results in the resulting '
    'visualization. Inputs are labeled with labels in the order that each '
    'is input. If there are fewer labels than inputs (or no labels), '
    'unnamed inputs are labeled numerically in sequential order. Extra '
    'labels are ignored.')

super_lca_desc = (
    '"super" finds the LCA consensus while giving preference to '
    'majority labels and collapsing substrings into superstrings. '
    'For example, when a more specific taxonomy does not '
    'contradict a less specific taxonomy, the more specific is '
    'chosen. That is, "g__Faecalibacterium; s__prausnitzii", '
    'will be preferred over "g__Faecalibacterium; s__"')

plugin.pipelines.register_function(
    function=evaluate_fit_classifier,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'reads_per_batch': _classify_parameters['reads_per_batch'],
        'n_jobs': _classify_parameters['n_jobs'],
        'confidence': _classify_parameters['confidence']},
    outputs=[('classifier', TaxonomicClassifier),
             ('evaluation', Visualization),
             ('observed_taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'sequences': 'Reference sequences to use for classifier '
                     'training/testing.',
        'taxonomy': 'Reference taxonomy to use for classifier '
                    'training/testing.'},
    parameter_descriptions={
        'reads_per_batch': _parameter_descriptions['reads_per_batch'],
        'n_jobs': _parameter_descriptions['n_jobs'],
        'confidence': _parameter_descriptions['confidence']},
    output_descriptions={
        'classifier': 'Trained naive Bayes taxonomic classifier.',
        'evaluation': 'Visualization of classification accuracy results.',
        'observed_taxonomy': 'Observed taxonomic label for each input '
                             'sequence, predicted by the trained classifier.'},
    name=('Evaluate and train naive Bayes classifier on reference sequences.'),
    description=(
        'Train a naive Bayes classifier on a set of reference sequences, then '
        'test performance accuracy on this same set of sequences. This '
        'results in a "perfect" classifier that "knows" the correct identity '
        'of each input sequence. Such a leaky classifier indicates the upper '
        'limit of classification accuracy based on sequence information '
        'alone, as misclassifications are an indication of unresolvable kmer '
        'profiles. This test simulates the case where all query sequences '
        'are present in a fully comprehensive reference database. To simulate '
        'more realistic conditions, see `evaluate_cross_validate`. THE '
        'CLASSIFIER OUTPUT BY THIS PIPELINE IS PRODUCTION-READY and can be '
        're-used for classification of other sequences (provided the '
        'reference data are viable), hence THIS PIPELINE IS USEFUL FOR '
        'TRAINING FEATURE CLASSIFIERS AND THEN EVALUATING THEM ON-THE-FLY.'),
    citations=[citations['bokulich2018optimizing']],
    migrated={'to_plugin': 'q2-feature-classifier',
              'from_distro': 'amplicon',
              'to_distro': 'amplicon', 'epoch': '2026.1'},
)


plugin.pipelines.register_function(
    function=evaluate_cross_validate,
    inputs={'sequences': FeatureData[Sequence],
            'taxonomy': FeatureData[Taxonomy]},
    parameters={
        'k': Int % Range(2, None),
        'random_state': Int % Range(0, None),
        'reads_per_batch': _classify_parameters['reads_per_batch'],
        'n_jobs': _classify_parameters['n_jobs'],
        'confidence': _classify_parameters['confidence']},
    outputs=[('expected_taxonomy', FeatureData[Taxonomy]),
             ('observed_taxonomy', FeatureData[Taxonomy]),
             ('evaluation', Visualization)],
    input_descriptions={
        'sequences': 'Reference sequences to use for classifier '
                     'training/testing.',
        'taxonomy': 'Reference taxonomy to use for classifier '
                    'training/testing.'},
    parameter_descriptions={
        'k': 'Number of stratified folds.',
        'random_state': 'Seed used by the random number generator.',
        'reads_per_batch': _parameter_descriptions['reads_per_batch'],
        'n_jobs': _parameter_descriptions['n_jobs'],
        'confidence': _parameter_descriptions['confidence']},
    output_descriptions={
        'expected_taxonomy': 'Expected taxonomic label for each input '
                             'sequence. Taxonomic labels may be truncated due '
                             'to k-fold CV and stratification.',
        'observed_taxonomy': 'Observed taxonomic label for each input '
                             'sequence, predicted by cross-validation.',
        'evaluation': 'Visualization of cross-validated accuracy results.'},
    name=('Evaluate DNA sequence reference database via cross-validated '
          'taxonomic classification.'),
    description=(
        'Evaluate DNA sequence reference database via cross-validated '
        'taxonomic classification. Unique taxonomic labels are truncated to '
        'enable appropriate label stratification. See the cited reference '
        '(Bokulich et al. 2018) for more details.'),
    citations=[citations['bokulich2018optimizing']],
    migrated={'to_plugin': 'q2-feature-classifier',
              'from_distro': 'amplicon',
              'to_distro': 'amplicon', 'epoch': '2026.1'},
)


palettes = ['Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired',
            'Accent', 'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c',
            'viridis', 'plasma', 'inferno', 'magma', 'cividis', 'terrain',
            'rainbow', 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
