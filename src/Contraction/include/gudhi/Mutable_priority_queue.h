/*
 * Mutable_priority_queue.h
 *
 *  Created on: Jun 3, 2014
 *      Author: dsalinas
 */

#ifndef MUTABLE_PRIORITY_QUEUE_H_
#define MUTABLE_PRIORITY_QUEUE_H_

namespace Gudhi{

template<typename EdgeHandle,typename CompareCost, typename ID> class Mutable_queue{

	boost::optional<Edge_handle> e = 	heap_PQ_->extract_top()
	data.set_PQ_handle(heap_PQ_->erase(edge,data.PQ_handle()));
	data.set_PQ_handle(heap_PQ_->push(edge));
	data.set_PQ_handle(heap_PQ_->update(edge,data.PQ_handle())) ;


//	get_data(*edge).reset_PQ_handle();
	heap_PQ_.reset( new PQ (size, Compare_cost(this), Undirected_edge_id(this) ) ) ;

	heap_PQ_->empty()
					auto edge = heap_PQ_->top();
};

}  // namespace GUDHI

#endif /* MUTABLE_PRIORITY_QUEUE_H_ */
